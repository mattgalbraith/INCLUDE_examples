############################################################
# linear_regression_example.R                              #
# Script original author: Matthew Galbraith                #
# version: 0.3  Date: 07_15_2022                           #
############################################################

if (!require("tidyverse")) install.packages("tidyverse"); library("tidyverse")
if (!require("tictoc")) install.packages("tictoc"); library("tictoc") # for timing
if (!require("rstatix")) install.packages("rstatix"); library("rstatix") # various function
if (!require("broom")) install.packages("broom"); library("broom") # for extracting model results etc
if (!require("ggforce")) install.packages("ggforce"); library("ggforce") # required for sina plots
if (!require("ggrepel")) install.packages("ggrepel"); library("ggrepel") # labels
if (!require("plotly")) install.packages("plotly"); library("plotly") # required for interactive plots
if (!require("here")) install.packages("here"); library("here") # managing dir paths

standard_colors <- c("Control" = "gray60", "T21" = "#009b4e")


# get list of files
metab_data_files <- list.files(path = "/sbgenomics/project-files/", pattern = "LCMS_Metabolomics.tsv.gz", full.names = TRUE)

# Concatenate individual sample-level files
tic()
metab_data <- msd_data_files %>% 
  map_dfr(~read_tsv(., id = "file")) %>% 
  mutate(filename = basename(file)) %>% 
  select(LabID, everything())
toc() # ~18.557 sec elapsed


# SHOULD BE IMPORTING CLINICAL DATA (Karyotype, Age, Sex, etc) AND JOINING HERE
# but we will cheat and get Karyotype (T21 status) information here from the "LabID"


# Prepare data for linear regression -----
# assumes  "metab_data" is in long format already joined with clinical data
# "LabID" is unique sample identifier; "Analyte" denotes each feature of interest; "Value" is the actual measurement; here "Karyotype" contains 2 levels: " and "T21"
regressions_dat <- metab_data %>%
  # select(LabID, Karyotype, Sex, Age, Analyte, Value) %>% 
  select(LabID, Analyte, Value) %>% 
  mutate(Karyotype = if_else(str_detect(LabID, "A"), "T21", "Control")) %>% # NEED TO REPLACE WITH ACTUAL CLINICAL DATA JOIN
  mutate(
    Karyotype = fct_relevel(Karyotype, c("Control", "T21")), # ensure factor levels in correct order
    #Sex = fct_relevel(Sex, c("Female", "Male")), # ensure factor levels in correct order
  ) %>% 
  group_by(Analyte, Karyotype) %>% # using both groupings here for categorical testing
  mutate(extreme = rstatix::is_extreme(log2(Value))) %>%
  ungroup() %>% 
  filter(extreme != TRUE) %>% # remove extreme outliers
  # I would usually also check here for:
  # 1) a minimum number of samples per group (eg 5) and 
  # 2) that there are >1 levels for categorical veriables of interest (prevents errors in regression step)
  nest(-Analyte) # nesting allows for easy testing of all features ~ at once


# Run simple linear regression for each feature with log2(Value) as outcome and "Karyotype" as predictor -----
tic("Running linear regressions for simple model...")
regressions_simple <- regressions_dat %>% 
  mutate(
    fit = map(data, ~ lm(log2(Value) ~ Karyotype, data = .x)),
    tidied = map(fit, broom::tidy), # see ?tidy.lm
    # glanced = map(fit, broom::glance), # see ?glance.lm # NOT NEEDED FOR DEMO
    # augmented = map(fit, broom::augment), # see ?augment.lm # NOT NEEDED FOR DEMO
  )
toc()
# Run linear regression for each feature with log2(Value) as outcome, "Karyotype" as predictor, and Age + Sex and nuisance variables -----
# THIS WILL CURRENTLY FAIL AS WE DO NOT HAVE AGE OR SEX VARIABLES!!!!
tic("Running linear regressions for multi model with Age + Sex...")
regressions_multi_SexAge <- regressions_dat %>% 
  mutate(
    fit = map(data, ~ lm(log2(Value) ~ Karyotype + Age + Sex, data = .x)),
    tidied = map(fit, broom::tidy), # see ?tidy.lm
    # glanced = map(fit, broom::glance), # see ?glance.lm # NOT NEEDED FOR DEMO
    # augmented = map(fit, broom::augment), # see ?augment.lm # NOT NEEDED FOR DEMO
    # vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term")) # NOT NEEDED FOR DEMO
  )
toc()


# Extract and plot results for simple lm: log2(Value) ~ Karyotype -----
lm_results_simple <- regressions_simple %>% 
  unnest(tidied) %>% 
  select(Analyte, term, estimate, p.value) %>% 
  group_by(Analyte) %>% 
  dplyr::summarize(
    Analyte = first(Analyte),
    log2_denom = first(estimate), # check for transformation and adjust accordingly
    log2_num = nth(estimate, n = 2) + log2_denom, # check for transformation and adjust accordingly
    log2FoldChange = nth(estimate, n = 2), # check for transformation and adjust accordingly; equivalent to difference between level 2 and level 1 ie y = B0 + B1x
    FoldChange = 2^log2FoldChange, # check for transformation and adjust accordingly
    pval = nth(p.value, n = 2)
  ) %>% 
  ungroup() %>% 
  arrange(pval) %>% 
  mutate(BHadj_pval = p.adjust(pval, method = "BH", n = length(pval))) %>% # correct for multiple hypothesis testing
  select(
    Analyte,
    log2_denom, # These should be relabelled according to factor of interest, eg "Control_mean"
    log2_num,  # These should be relabelled according to factor of interest, eg "T21_mean"
    FoldChange,
    log2FoldChange,
    pval,
    BHadj_pval,
    everything()
  )
#
# Extract and plot results for multi lm: log2(Value) ~ Karyotype + Age + Sex -----
lm_results_multi_SexAge <- regressions_multi_SexAge %>% 
  unnest(tidied) %>% 
  select(Analyte, term, estimate, p.value) %>% 
  group_by(Analyte) %>% 
  dplyr::summarize(
    Analyte = first(Analyte),
    log2_denom = first(estimate), # check for transformation and adjust accordingly
    log2_num = nth(estimate, n = 2) + log2_denom, # check for transformation and adjust accordingly
    log2FoldChange = nth(estimate, n = 2), # check for transformation and adjust accordingly; equivalent to difference between level 2 and level 1 ie y = B0 + B1x
    FoldChange = 2^log2FoldChange, # check for transformation and adjust accordingly
    pval = nth(p.value, n = 2)
  ) %>% 
  ungroup() %>% 
  arrange(pval) %>% 
  mutate(BHadj_pval = p.adjust(pval, method = "BH", n = length(pval))) %>% 
  select(
    Analyte,
    log2_denom, # These should be relabelled according to factor of interest, eg "Control_mean"
    log2_num,  # These should be relabelled according to factor of interest, eg "T21_mean"
    FoldChange,
    log2FoldChange,
    pval,
    BHadj_pval,
    everything()
  )
#


# function for Volcano plot of lm results
volcano_plot_lab_lm <- function(res, title = "", 
                                subtitle = "",
                                y_lim = c(0, NA)){
  theme_set(theme_gray(base_size=12, base_family="Arial") +
              theme(panel.border=element_rect(colour="black", fill="transparent"), 
                    plot.title=element_text(face="bold", hjust=0), # lineheight=.8, size=20,
                    axis.text=element_text(color="black", size=14), 
                    axis.text.x=element_text(angle=0, hjust=0.5),
                    panel.background=element_blank(),
                    panel.grid=element_blank(),
                    plot.background=element_blank(),
                    strip.background = element_blank(), # facet label borders
                    legend.key=element_blank(), legend.background=element_blank() # remove grey bg from legend
              )
  )
  res <- res %>% 
    mutate(
      color = if_else(BHadj_pval < 0.1, "padj < 0.1", "All")
    )
  # get max for x-axis
  x_lim <- res %>% 
    summarize(max = max(log2(FoldChange), na.rm = TRUE), min = min(log2(FoldChange), na.rm = TRUE)) %>% 
    abs() %>% 
    max() %>% 
    ceiling()
  # generate plot
  res %>% 
    ggplot(aes(log2(FoldChange), -log10(BHadj_pval), color = color, label = Analyte)) + # need label = "Analyte" for reature names to work in ggplotly
    geom_hline(yintercept = -log10(0.1), linetype = 2) + 
    geom_vline(xintercept = 0, linetype = 2) + 
    geom_point() + 
    scale_color_manual(values = c("padj < 0.1" = "red", "All" = "black")) + 
    xlim(-x_lim, x_lim) + # set x-axis to be symmetrical
    ylim(y_lim) +
    # labels for top features by fold-change and/or adjusted p
    geom_text_repel(data = res %>% filter(!is.na(BHadj_pval)) %>% slice_max(order_by = FoldChange, n = 3), aes(label = Analyte), min.segment.length = 0, show.legend = FALSE, nudge_x = 1, nudge_y = 0.1) +
    geom_text_repel(data = res %>% filter(!is.na(BHadj_pval)) %>% slice_min(order_by = FoldChange, n = 3), aes(label = Analyte), min.segment.length = 0, show.legend = FALSE, nudge_x = -1, nudge_y = 0.1) +
    # geom_text_repel(data = res %>% filter(!is.na(BHadj_pval)) %>% slice_min(order_by = BHadj_pval, n = 3), aes(label = Analyte), min.segment.length = 0, show.legend = FALSE, nudge_x = -0, nudge_y = 0.1) +
    theme(aspect.ratio=1.2) +
    labs(
      title = title,
      subtitle = subtitle
    )
} # end of function

# Volcano plot for simple model ----
lm_results_simple %>% 
  volcano_plot_lab_lm(
    title="Differential abundance in T21 vs. Controls",
    # subtitle with number significant up/down:
    subtitle = paste0("[Down: ",(.) %>% filter(BHadj_pval < 0.1 & FoldChange <1) %>% nrow(), "; Up: ",(.) %>% filter(BHadj_pval < 0.1 & FoldChange >1) %>% nrow(), "]")
  )
ggsave("/sbgenomics/output-files/volcano.png")
# can pass to plotly to make interactive plot
plotly::ggplotly()


# Modified Sina + Boxplot to display for individual features ----
# No adjustment for nuisance variables
regressions_dat %>% 
  unnest(data) %>% # ensures using exact same data as model
  filter(Analyte %in% c("Taurolithocholic acid", "5(S)-HETE")) %>% # FILTER TO FEATURE(S) OF INTEREST
  mutate(Analyte = fct_relevel(Analyte, c("Taurolithocholic acid", "5(S)-HETE"))) %>% # if >1 feature relevel to control plotting order
  ggplot(aes(Karyotype, log2(Value), color = Karyotype)) + # SET X AND COLOR TO VARIABLE OF INTEREST
  geom_sina() + # HORIZONTAL JITTERING OF POINTS BASED ON LOCAL DENSITY
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) + # WHISKERS AND OUTLIERS TURNED OFF BECAUSE SINA SHOWS ALL DATA POINTS
  facet_wrap(~ Analyte, scales = "free_y") + # EVEN FOR SINGLE FEATURE ENSURES NICE FEATURE LABEL ABOVE PLOT
  scale_color_manual(values = standard_colors) + # DEFINED AT TOP
  labs(
    title = "Plasma metabolites: T21 vs Control",
  )
ggsave("/sbgenomics/output-files/sina.png")




