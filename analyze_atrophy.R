library(dplyr)
library(stringr)
library(ggpubr)
library(qqplotr)
library(reshape2)
library(grid)
library(cowplot)


RESULTS_DIR <- './'
STATS_DIR <- 'DL+DiReCT_results_stats2table'


theme_set(theme_classic())
scale_colour_discrete_orig <- scale_colour_discrete
scale_colour_discrete <- function(...) { scale_colour_brewer(..., palette="Dark2") }

SHAPE_LOOKUP <- setNames(c(19, 3),
                         c('1.5', '3'))

# Read stats2table like files from FreeSurfer or DL+DiReCT
read_stats <- function(dataset, stat_files, sep='\t', base_dir=RESULTS_DIR) {
  all_stats <- data.frame()
  for(idx in c(1:length(stat_files))) {
    filename <- sprintf('%s/%s/%s', base_dir, dataset, stat_files[idx])
    print(filename)
    stats <- read.csv(filename, header = TRUE, sep=sep, stringsAsFactors=FALSE)
    colnames(stats)[1] <- 'SUBJECT_ID'
    
    if(nrow(all_stats) == 0) {
      all_stats <- stats
    } else {
      all_stats <- merge(all_stats, stats, by='SUBJECT_ID') %>% dplyr::select(-ends_with('.x')) %>% select(-ends_with('.y'))
    }
  }
  
  return(all_stats)
}

# read DICOM metadata
metadata_msc <- read.csv(paste(RESULTS_DIR, 'metadata_dicom.csv', sep='/'), header = TRUE, sep = ";", stringsAsFactors=TRUE)
metadata_msc$SUBJECT_ID <- as.character(metadata_msc$SUBJECT_ID)
metadata_msc$SOURCE_SUBJECT <- gsub('([^0-9]+[0-9]+)_.*', '\\1', metadata_msc$SUBJECT_ID)

# Assuming contrast-enhanced MRI end with _CE
metadata_msc$IMAGE <- ifelse(endsWith(metadata_msc$SUBJECT_ID, '_CE'), 'CE', 'Native')

# e.g. treat TI=1.17 as TI=1.1
metadata_msc$InversionTimeRounded <- floor(metadata_msc$InversionTime*10)/10
metadata_msc$InversionTime <- factor(metadata_msc$InversionTime)
metadata_msc$MagneticFieldStrength <- factor(round(metadata_msc$MagneticFieldStrength, 1))

# read mapping and clinical data (with EventDate)
mapping <- read.csv(paste(RESULTS_DIR, 'mapping.csv', sep='/'), header = TRUE, sep = ";", stringsAsFactors=TRUE)
clinical_data <- read.csv(paste(RESULTS_DIR, 'clinical_data.csv', sep='/'), header = TRUE, sep = ",", stringsAsFactors=TRUE)

meta_add <- merge(mapping, clinical_data, by='PatientID', all.x=TRUE)
meta_add$StudyDate <- as.Date(as.character(meta_add$StudyDate), format='%Y%m%d')
meta_add$EventDate <- as.Date(as.character(meta_add$EventDate), format='%Y%m%d')
meta_add$PatientBirthDate <- as.Date(as.character(meta_add$PatientBirthDate), format='%Y%m%d')

metadata_msc <- merge(metadata_msc, meta_add %>% select(SUBJECT_ID, StudyDate, EventDate, PatientBirthDate), by='SUBJECT_ID', all.x=TRUE)

metadata_msc$AgeAtTest <- as.numeric(round((metadata_msc$EventDate - metadata_msc$PatientBirthDate)/365.25, 2))
metadata_msc$AgeToday <- as.numeric(round((Sys.Date() - metadata_msc$PatientBirthDate)/365.25, 2))
metadata_msc$TimeRelativeToEvent <- as.numeric(round(metadata_msc$AGE - metadata_msc$AgeAtTest, 2))
metadata_msc$EventPeriod <- factor(ifelse(metadata_msc$AGE <= metadata_msc$AgeAtTest, 'Pre', 'Post'), levels=c('Pre', 'Post'))

EXCLUDES = c()
metadata_msc <- metadata_msc %>% filter(!SUBJECT_ID %in% EXCLUDES)


# identify matching sequences
metadata_msc$IS_SEQ_MATCH = FALSE
metadata_msc$INCLUDE = FALSE
for(subj in unique(metadata_msc$SOURCE_SUBJECT)) {
  for(img_type in unique(metadata_msc$IMAGE)) {
    df_subj <- metadata_msc %>% filter(SOURCE_SUBJECT == subj & IMAGE == img_type)
    
    # post-event, TI=1.1
    candidates_post <- df_subj %>% filter(EventPeriod == 'Post' & InversionTimeRounded == 1.1)
    
    params <- candidates_post %>% group_by(MagneticFieldStrength, RepetitionTime, InversionTimeRounded) %>% select(MagneticFieldStrength, RepetitionTime, InversionTimeRounded) %>% mutate(N=n()) %>% unique() %>% arrange(desc(N))
    if(nrow(params) > 0) {
      match_params <- params[1, ]
      
      print(sprintf('%s/%s: %s / %s / %s', subj, img_type, match_params$MagneticFieldStrength, match_params$RepetitionTime, match_params$InversionTimeRounded))
      
      metadata_msc[metadata_msc$SOURCE_SUBJECT == subj & metadata_msc$IMAGE == img_type & 
                   metadata_msc$MagneticFieldStrength == match_params$MagneticFieldStrength &
                   metadata_msc$InversionTimeRounded == match_params$InversionTimeRounded &
                   abs(metadata_msc$RepetitionTime - match_params$RepetitionTime) <= 0.11, 'IS_SEQ_MATCH'] <- TRUE
      
      # check if at least 3 matching pre-event sequences
      if(nrow(metadata_msc %>% filter(SOURCE_SUBJECT == subj & IMAGE == img_type & IS_SEQ_MATCH == TRUE & EventPeriod == 'Pre')) >= 3) {
        metadata_msc[metadata_msc$SOURCE_SUBJECT == subj & metadata_msc$IMAGE == img_type & metadata_msc$IS_SEQ_MATCH == TRUE, 'INCLUDE'] <- TRUE
      }
    } else {
      print(sprintf('%s/%s: NA', subj, img_type))
    }
  }
}

highlight_subj <- function(p, img_type) {
  subj_includes <- c(metadata_msc %>% filter(IMAGE == img_type & INCLUDE == TRUE) %>% select(SOURCE_SUBJECT) %>% unique())$SOURCE_SUBJECT
  g <- ggplot_gtable(ggplot_build(p))
  stripr <- which(grepl('strip-', g$layout$name))
  for (i in stripr) {
    lbl <- g$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$label
    if(!is.null(lbl) && lbl %in% subj_includes) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- 'lightblue'
    }
  }
  return(ggdraw()+draw_plot(g))
}

lim <- c(min(metadata_msc$RepetitionTime)*0.9, max(metadata_msc$RepetitionTime)*1.1)
labels <- metadata_msc %>% group_by(SOURCE_SUBJECT) %>% mutate(AgeAtTest=AgeAtTest, MaxAge=AgeToday) %>% select(SOURCE_SUBJECT, AgeAtTest, MaxAge) %>% unique()
plot_subj_params <- function(df, type) {
  p <- ggplot(df %>% filter(IMAGE == type), aes(x=AGE, y=RepetitionTime, color=InversionTime, shape=MagneticFieldStrength, size=IS_SEQ_MATCH))+
    geom_jitter(height=0, alpha=0.5)+
    scale_size_manual(values=c(1, 2.5))+
    scale_shape_manual(values=SHAPE_LOOKUP[levels(df$MagneticFieldStrength)], labels = levels(df$MagneticFieldStrength))+
    geom_rect(data=labels, aes(xmin = AgeAtTest, xmax = MaxAge, ymin = -Inf, ymax = Inf, alpha=0.1), inherit.aes = FALSE, fill='lightgrey') +
    scale_alpha(guide = 'none')+labs(shape='Tesla', color='TI', size='Matching Seq.')+#theme(legend.justification=c(1,0), legend.position=c(1,0))+
    facet_wrap(~SOURCE_SUBJECT, scales='free_x')+ylim(lim)#+labs(title=type)
  
  return(highlight_subj(p, type))
}

p1 <- plot_subj_params(metadata_msc, 'CE')
p2 <- plot_subj_params(metadata_msc, 'Native')
ggarrange(plotlist=list(p1, p2), nrow=1, ncol=2)
#ggsave(filename = 'subj_params.pdf', width=16, height = 10, device=cairo_pdf)
#ggsave(p2, filename = 'subj_params_native.pdf', width=10, height = 10, device=cairo_pdf)


df_results <- read_stats(STATS_DIR, c('aseg_stats_volume.txt', 'lh.aparc_stats_thickness.txt', 'rh.aparc_stats_thickness.txt', 'lh.aparc_stats_volume.txt', 'rh.aparc_stats_volume.txt'), sep='\t')
df_results$MeanThickness <- (df_results$lh_MeanThickness_thickness + df_results$rh_MeanThickness_thickness)/2
df_results <- merge(df_results, metadata_msc, by='SUBJECT_ID')

df_results$CortexVol <- df_results$TotalCortexVol
df_results$bila_parahippocampal_volume <- df_results$lh_parahippocampal_volume + df_results$rh_parahippocampal_volume


plot_subj <- function(df, metric, matches_only=TRUE, multiplier=1, ytitle=NULL) {
  ytitle <- ifelse(is.null(ytitle), metric, ytitle)
  df_pre <- df %>% filter(EventPeriod == 'Pre')
  if(matches_only) {
    df_pre <- df_pre %>% filter(IS_SEQ_MATCH == TRUE & INCLUDE==TRUE)
  }
  p <- ggplot(df, aes(x=AGE, y=df[, metric]*multiplier, shape=MagneticFieldStrength, size=IS_SEQ_MATCH, color=INCLUDE))+
    geom_point()+
    scale_size_manual(values=c(1, 2.5))+
    scale_shape_manual(values=SHAPE_LOOKUP[levels(df$MagneticFieldStrength)], labels = levels(df$MagneticFieldStrength))+
    geom_rect(data=labels, aes(xmin = AgeAtTest, xmax = MaxAge, ymin = -Inf, ymax = Inf, alpha=0.1), inherit.aes = FALSE, fill='lightgrey')+
    geom_smooth(aes(x=AGE, y=df_pre[, metric]*multiplier), data=df_pre, inherit.aes = FALSE, method = 'lm', se=FALSE, linetype='dashed', size=0.4, fullrange=TRUE)+
    scale_alpha(guide = 'none')+labs(shape='Tesla', color='Include', size='Matching Seq.')+#theme(legend.justification=c(1,0), legend.position=c(1,0))+
    facet_wrap(~SOURCE_SUBJECT, scales='free_x')+ylab(ytitle)+xlab('Age [years]')#+ylim(lim)+labs(title=type)
  return(highlight_subj(p, unique(df$IMAGE)))
}

p1 <- plot_subj(df_results %>% filter(IMAGE == 'CE'), 'MeanThickness', matches_only = T)
p2 <- plot_subj(df_results %>% filter(IMAGE == 'Native'), 'MeanThickness', matches_only = T)
ggarrange(plotlist=list(p1, p2), nrow=1, ncol=2)


p1 <- plot_subj(df_results %>% filter(IMAGE == 'CE'), 'TotalCortexVol', matches_only = T, multiplier = 0.001, ytitle='Cortical GM Volume [ml]')
p2 <- plot_subj(df_results %>% filter(IMAGE == 'Native'), 'TotalCortexVol', matches_only = T, multiplier = 0.001, ytitle='Cortical GM Volume [ml]')
ggarrange(plotlist=list(p1, p2), nrow=1, ncol=2)
#ggsave(filename = 'vol.pdf', width=16, height = 10, device=cairo_pdf)
#ggsave(p2, filename = 'vol_native.pdf', width=10, height = 10, device=cairo_pdf)


p1 <- plot_subj(df_results %>% filter(IMAGE == 'Native'), 'lh_parahippocampal_volume', matches_only = T, multiplier = 0.001)
p2 <- plot_subj(df_results %>% filter(IMAGE == 'Native'), 'rh_parahippocampal_volume', matches_only = T, multiplier = 0.001)
p3 <- plot_subj(df_results %>% filter(IMAGE == 'Native'), 'bila_parahippocampal_volume', matches_only = T, multiplier = 0.001)
ggarrange(plotlist=list(p1, p2, p3), nrow=3, ncol=1)
#ggsave(filename = 'vol_parahippocampal_native.pdf', width=16, height = 30, device=cairo_pdf)


calc_residuals <- function(df, metric, include_only=TRUE) {
  df[, paste0(metric, '.residual')] <- NA
  df[, paste0(metric, '.sresidual')] <- NA
  
  for(img in unique(df$IMAGE)) {
    for(subj in unique(df$SOURCE_SUBJECT)) {
      df_subset <- df %>% filter(IMAGE == img & SOURCE_SUBJECT == subj)
      if(include_only) {
        df_subset <- df_subset %>% filter(INCLUDE == TRUE)
      }
      df_subset_pre <- df_subset %>% filter(EventPeriod == 'Pre')
      
      if(nrow(df_subset_pre) > 0) {
        subj.lm <- lm(as.formula(paste0(metric, '~AGE')), data=df_subset_pre)
        residuals <- df_subset[, metric] - predict.lm(subj.lm, df_subset)
        
        df[df$SUBJECT_ID %in% df_subset$SUBJECT_ID, paste0(metric, '.residual')] <- residuals
        df[df$SUBJECT_ID %in% df_subset$SUBJECT_ID, paste0(metric, '.sresidual')] <- (residuals - mean(residuals)) / sd(residuals)
      } else {
        #print(sprintf('Excluding: %s/%s', subj, img))
      }
    }
  }
  
  return(df)
}


rois <- names(df_results)[grepl('(lh|rh|bila)_.*(volume|thickness)$|^MeanThickness$|^CortexVol$', names(df_results))]
for(roi in rois) {
  df_results <- calc_residuals(df_results, roi)
}


boxplot_res <- function(df, metric, title=NULL) {
  metric <- paste0(metric, '.sresidual')
  lim <- c(min(df[, metric], na.rm=TRUE), max(df[, metric], na.rm=TRUE)*1.2)
  p <- ggplot(df %>% filter(INCLUDE == TRUE), aes_string(x='EventPeriod', y=metric, color='EventPeriod'))+
    geom_boxplot()+ #(outlier.colour = NA)+
    stat_compare_means(method='t.test', label.x.npc='center', vjust=-0.0, comparisons = list(c('Pre', 'Post')))+
    ylab('')+
    labs(title=' ')+ # Note: Need empty title to make sure plot is aligned with left scatter plot
    theme(legend.position='none')+ylim(lim)
  
  if(length(unique(df$IMAGE)) > 1) {
    p <- p+
      facet_wrap(~IMAGE, scales='free_x')
  }
  
  return(p)
}

plot_res <- function(df, metric, title=NULL) {
  metric <- paste0(metric, '.sresidual')
  lim <- c(min(df[, metric], na.rm=TRUE), max(df[, metric], na.rm=TRUE)*1.2)
  print(lim)
  p <- ggplot(df %>% filter(INCLUDE == TRUE), aes_string(x='TimeRelativeToEvent', y=metric, color='EventPeriod'))+
    geom_point()+
    scale_y_continuous(expand = expansion(mult = c(0, 0.3)))+
    geom_vline(xintercept = 0, size=0.3, linetype='dashed')+
    xlab('Years relative to event')+ylab('Stand. Residuals')+labs(title=title)+
    theme(legend.position='none', plot.title = element_text(size=12, hjust=0.5, face='bold'))+ylim(lim)
  
  if(length(unique(df$IMAGE)) > 1) {
    p <- p+
      facet_wrap(~IMAGE)#, scales='free_x')
  }
  
  return(p)
}


p1 <- plot_res(df_results %>% filter(IMAGE == 'Native'), 'CortexVol', 'Cortical GM Volume')
p2 <- boxplot_res(df_results %>% filter(IMAGE == 'Native'), 'CortexVol')
p3 <- plot_res(df_results %>% filter(IMAGE == 'Native'), 'bila_parahippocampal_volume', 'Parahippocampal Gyri Volume')
p4 <- boxplot_res(df_results %>% filter(IMAGE == 'Native'), 'bila_parahippocampal_volume')
ggarrange(plotlist=list(p1, p2, p3, p4), nrow=2, ncol=2, widths = c(1, 0.2))
#ggsave(filename = 'results_cortex.pdf', width=12, height = 6, device=cairo_pdf)


p1 <- plot_res(df_results %>% filter(IMAGE == 'Native'), 'lh_parahippocampal_volume', 'Left Parahippocampal Gyrus Volume')
p2 <- boxplot_res(df_results %>% filter(IMAGE == 'Native'), 'lh_parahippocampal_volume')
p3 <- plot_res(df_results %>% filter(IMAGE == 'Native'), 'rh_parahippocampal_volume', 'Right Parahippocampal Gyrus Volume')
p4 <- boxplot_res(df_results %>% filter(IMAGE == 'Native'), 'rh_parahippocampal_volume')
p5 <- plot_res(df_results %>% filter(IMAGE == 'Native'), 'bila_parahippocampal_volume', 'Bilateral Parahippocampal Gyri Volume')
p6 <- boxplot_res(df_results %>% filter(IMAGE == 'Native'), 'bila_parahippocampal_volume')
ggarrange(plotlist=list(p1, p2, p3, p4, p5, p6), nrow=3, ncol=2, widths = c(1, 0.2))
#ggsave(filename = 'results_parahippocampal.pdf', width=12, height = 9, device=cairo_pdf)

