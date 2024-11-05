rm(list = ls())
library(tidyverse)
library(latex2exp)
library(scales)
library(xtable)

# set my ggplot theme
mytheme <- theme_minimal(base_size = 19) + 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(color = "grey30"),
        axis.text = element_text(color = "black"),
        axis.title.x = element_text(margin = margin(t = 12, r = 0, b = 0, l = 0)))

theme_set(mytheme)

overall_df_dml <- readRDS("01_dml_log_1000.RDS")
overall_df_rwr <- readRDS("02_rwr_log_1000.RDS")

table1_dml <-  overall_df_dml %>%
  filter(estimand == "tau0" | str_detect(estimand, "theta")) %>%
  mutate_at(vars(est, se), ~ formatC(.x, format = "f", digits = 3)) %>%
  transmute(estimand = fct_relevel(estimand, c("tau0", "theta0", "theta1", "theta2", "theta3")),
            estimate = paste0(est, " (", se, ")")) %>%
  arrange(estimand) %>%
  pivot_wider(names_from = estimand, values_from = estimate)  %>% mutate(" " = " DML") %>% dplyr::select(" ", everything())

table1_rwr <-  overall_df_rwr %>%
  filter(estimand == "tau0" | str_detect(estimand, "theta")) %>%
  mutate_at(vars(est, se), ~ formatC(.x, format = "f", digits = 3)) %>%
  # filter(estimand %in% c("ate", "cde0", "pie", "refint", "medint", "atm")) %>%
  transmute(estimand = fct_relevel(estimand, c("tau0", "theta0", "theta1", "theta2", "theta3")),
            estimate = paste0(est, " (", se, ")")) %>%
  arrange(estimand) %>%
  pivot_wider(names_from = estimand, values_from = estimate) %>% mutate(" " = " RWR") %>% dplyr::select(" ", everything())

table1 <- bind_rows(table1_dml, table1_rwr)


##################################################
#################### FIGURE 2 ####################
##################################################

thetas_out_df <- 
  bind_rows(overall_df_dml %>%
              filter(estimand == "tau0" | str_detect(estimand, "theta")) %>% mutate(estimator = "DML"),
            overall_df_rwr %>%
              filter(estimand == "tau0" | str_detect(estimand, "theta")) %>% mutate(estimator = "RWR")) %>%
  mutate(estimand = factor(estimand, 
                           levels = rev(c("tau0", "theta0", "theta1", "theta2", "theta3")),
                           # labels = unname(TeX(c("Gd. School \nAttendance ($\\Delta_1$)", "BA \nCompletion ($A_3$)",
                           #                       "College \nAttendance ($A_2$)", "HS \nCompletion ($A_1$)", "test")))),
                           labels = rev(c(expression(paste("Total Effect (", tau[`0`], ")")),
                                          expression(paste("Direct Effect (", Delta[`0`], ")")),
                                          expression(paste("via College Attendance (", theta[`1`], ")")),
                                          expression(paste("via College Completion (", theta[`2`], ")")),
                                          expression(paste("via Grad. School Attendance (", theta[`3`], ")"))))))

ggplot(thetas_out_df, aes(x = estimand, y = est, color = estimator, shape = estimator)) +
  geom_pointrange(aes(ymin = est - 1.96 * se,  ymax = est + 1.96 * se),
                  position = position_dodge(width = - 0.5), size = .9) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual("", values = c("blueviolet", "red"), labels = unname(TeX(c("DML", "RWR"))))  +  
  scale_shape("", labels = unname(TeX(c("DML", "RWR"))))  + 
  scale_color_manual("", values = c("blueviolet", "red"), labels = unname(TeX(c("DML", "RWR"))))  +  
  scale_x_discrete("", labels = parse_format()) +
  scale_y_continuous("Effect of HS Completion on Log Earnings") +
  coord_flip()


ggsave("Output/fig2.png", width = 9, height = 6)

##################################################
################### TABLE 1 ######################
##################################################

table2a_dml <-  overall_df_dml %>%
  filter(!str_detect(estimand, "tau") , !str_detect(estimand, "theta")) %>%
  mutate_at(vars(est, se), ~ formatC(.x, format = "f", digits = 3)) 

table2b_dml <- table2a_dml %>%
  transmute(estimand = fct_relevel(estimand, c("delta0", "delta1", "delta2", "delta3", "pi1", "pi2", "pi3", "eta1", "eta2", "eta3")),
            estimate = paste0(est)) %>%
  arrange(estimand) %>%
  pivot_wider(names_from = estimand, values_from = estimate)

table2_dml <- table2b_dml %>% rbind(paste0("(", table2a_dml$se, ")"))%>% mutate(" " = c(" DML", "")) %>% select(" ", everything())

table2a_rwr <-  overall_df_rwr %>%
  filter(!str_detect(estimand, "tau") , !str_detect(estimand, "theta")) %>%
  mutate_at(vars(est, se), ~ formatC(.x, format = "f", digits = 3)) 

table2b_rwr <- table2a_rwr %>%
  transmute(estimand = fct_relevel(estimand, c("delta0", "delta1", "delta2", "delta3", "pi1", "pi2", "pi3", "eta1", "eta2", "eta3")),
            estimate = paste0(est)) %>%
  arrange(estimand) %>%
  pivot_wider(names_from = estimand, values_from = estimate)

table2_rwr <- table2b_rwr %>% rbind(paste0("(", table2a_rwr$se, ")"))%>% mutate(" " = c(" RWR", "")) %>% select(" ", everything())

table2 <- bind_rows(table2_dml, table2_rwr)


table2.x <- xtable(table2, type = "latex")

# Add to row command
addtorow <- list()
addtorow$pos <- list(0,0,0,4)
hline <- paste0("\\hline")
line1 <- paste0("& $\\Delta_0$ & $\\Delta_1$ & $\\Delta_2$ & $\\Delta_3$ & $\\pi_1$ & $\\pi_2$ & $\\pi_3$ & $\\eta_1 & $\\eta_2$ & $\\eta_3$  \\\\")
hline <- paste0("\\hline")
#line2 <- paste0("\\multicolumn{11}{l}{\\small  Note: Numbers in parentheses are estimates of standard errors, adjusted for multiple imputation}", '\\\\')
#line3 <- paste0("\\multicolumn{11}{l}{\\small  via Rubin`s (1987) method.}", '\\\\')
addtorow$command <- c(hline, line1, hline, hline)

print(table2.x, add.to.row=addtorow, include.colnames=F,  include.rownames = F, hline.after = c())
print(table2.x, add.to.row=addtorow, include.colnames=F,  include.rownames = F, size="\\fontsize{11pt}{14pt}\\selectfont", hline.after = c(),
      file = "Output/table1.tex")

