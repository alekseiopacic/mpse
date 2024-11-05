library(haven)
library(splines)
library(survey)
library(rbw)
library(gbm)
library(ranger)
library(rlang)
library(glmnet)
library(kernlab)
library(nnet)
library(rsample)
library(caret)
library(CBPS)
library(quantreg)
library(KernSmooth)
library(mgcv)
library(nprobust)
library(mice)
library(SuperLearner)
library(tidyverse)
library(gridExtra)
# library(rootSolve)
library(Hmisc)

source("zmisc.R")
set.seed(02138)
# load data ; construct treatment and mediator variables



#load("Output/04_chdearn_overall_dml.RData")

var <- "chdearn2" # , "chdearn_old")  # 32 to 36, 30-33
nlsy97_full <-  readRDS(file = paste0("Samples/", var, "drop_full.RDS")) %>%
  map( ~ mutate(.x, 
                hs = if_else(educAB2_grad >= 2, 1, 0),
                college = if_else(educAB2_grad >= 4, 1, 0),
                ba = if_else(educAB2_grad >= 5, 1, 0),
                gschool = if_else(educAB2_grad >= 6, 1, 0)))  %>%
  map( ~ mutate(.x, earnings = log(!!sym(var) + 1000),
                afqt3_bi = if_else(afqt3 >= median(afqt3, na.rm = T), 1, 0)))



overall_df <- readRDS("01_dml_log_1000.RDS") %>%
  filter(str_detect(estimand, 'delta|tau|pi'))
mytheme <- theme_minimal(base_size = 14) + 
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 16.5),
        plot.caption = element_text(color = "grey30"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
theme_set(mytheme)

########## Get the estimate of interventional gap for AA-equalization


#####The curve that renders interventional disparity estimate 0: 0.026 = 0.18 \gamma\delta
#####The curve that renders interventional disparity estimate insignificant: 0.013 = 0.18 \gamma\delta

###########################################
############ Observed conditional associations ############
#############################################
a <- expr(hs)
m1 <- expr(college)
m2 <- expr(ba)
m3 <- expr(gschool)
y <- expr(earnings)
w <- expr(weight)
id <- expr(id)

x <- nlsy97_full[[1]] %>%
  dplyr::select(female, black, hispan, parinc, paredu, parast, livebio, father, # no other! - i filtered on ths
                afqt3, gpa, children18, substind, delinqind, rural97, south97,
                frcoll_75, frcoll_90, schstolen, schthreat, schfight, afqt3_bi) %>%
  names()

z <- nlsy97_full[[1]] %>%
  dplyr::select( stem, college_gpa) %>%
  names()

list_t0 <- list()
for (var in x) {
  
  x_new <- eval(parse_expr(paste0(" nlsy97_full[[1]] %>% dplyr::select(", paste(x, collapse = ", "), ") %>% dplyr::select(-", var, ")"))) %>% names()
  yx_mod <- svyglm(as.formula(paste(y, " ~ ", paste(c(x, a), collapse= "+"))), design = svydesign(ids = ~1, data = nlsy97_full[[1]], weights = ~weight))
  ax_mod <- svyglm(as.formula(paste(var, " ~ ", paste(c(x_new, a), collapse= "+"))), design = svydesign(ids = ~1, data = nlsy97_full[[1]], weights = ~weight))
  
  list_t0[[""]] <- tibble(name = var,
                          ymod = coef(yx_mod)[[var]],
                          amod = coef(ax_mod)[["hs"]])
}

list_t1d0 <- list()
for (var in x) {
  
  newd <- nlsy97_full[[1]] %>% filter(hs == 1)
  x_new <- eval(parse_expr(paste0(" newd %>% dplyr::select(", paste(x, collapse = ", "), ") %>% dplyr::select(-", var, ")"))) %>% names()
  yx_mod <- svyglm(as.formula(paste(y, " ~ ", paste(c(x, m1), collapse= "+"))), design = svydesign(ids = ~1, data = newd, weights = ~weight))
  ax_mod <- svyglm(as.formula(paste(var, " ~ ", paste(c(x_new, m1), collapse= "+"))), design = svydesign(ids = ~1, data = newd, weights = ~weight))
  
  list_t1d0[[""]] <- tibble(name = var,
                            ymod = coef(yx_mod)[[var]],
                            amod = coef(ax_mod)[["college"]])
}


list_t2d1 <- list()
for (var in c(x,z)) {
  
  newd <- nlsy97_full[[1]] %>% filter(college == 1)
  x_new <- eval(parse_expr(paste0(" newd %>% dplyr::select(", paste(c(x,z), collapse = ", "), ") %>% dplyr::select(-", var, ")"))) %>% names()
  yx_mod <- svyglm(as.formula(paste(y, " ~ ", paste(c(x, z, m2), collapse= "+"))), design = svydesign(ids = ~1, data = newd, weights = ~weight))
  ax_mod <- svyglm(as.formula(paste(var, " ~ ", paste(c(x_new, m2), collapse= "+"))), design = svydesign(ids = ~1, data = newd, weights = ~weight))
  
  list_t2d1[[""]] <- tibble(name = var,
                            ymod = coef(yx_mod)[[var]],
                            amod = coef(ax_mod)[["ba"]])
}

list_t3d2 <- list()
for (var in c(x,z)) {
  
  newd <- nlsy97_full[[1]] %>% filter(ba == 1)
  x_new <- eval(parse_expr(paste0(" newd %>% dplyr::select(", paste(c(x,z), collapse = ", "), ") %>% dplyr::select(-", var, ")"))) %>% names()
  yx_mod <- svyglm(as.formula(paste(y, " ~ ", paste(c(x, z, m3), collapse= "+"))), design = svydesign(ids = ~1, data = newd, weights = ~weight))
  ax_mod <- svyglm(as.formula(paste(var, " ~ ", paste(c(x_new, m3), collapse= "+"))), design = svydesign(ids = ~1, data = newd, weights = ~weight))
  
  list_t3d2[[""]] <- tibble(name = var,
                            ymod = coef(yx_mod)[[var]],
                            amod = coef(ax_mod)[["gschool"]])
}


out = lapply(list(list_t0, list_t1d0, list_t2d1, list_t3d2), function(t) { do.call(rbind,t )}) %>% do.call(rbind,. ) %>% mutate(prod = ymod * amod)
print(out,n=110)

out = lapply(lapply(list(list_t0, list_t1d0, list_t2d1, list_t3d2), function(t) { do.call(rbind,t )}), function(z) { z %>% mutate(prod = ymod * amod)})


# #Gender, gamma (E[Y=1|X=1]-E[Y=1|X=0]) and delta (P[X=1|A=1]-P[X=1|A=0])
# # yxgender <- svyglm(comp ~ female, design = svydesign(ids = ~1, data = nlsy97_full[[1]], weights = ~weight))
# # axgender <- svyglm(female ~ sel, design = svydesign(ids = ~1, data = nlsy97_full[[1]], weights = ~weight))
# 
# #father absence, gamma and delta
# yxfather <- svyglm(comp ~ father, design = svydesign(ids = ~1, data = nlsy97_full[[1]], weights = ~weight))
# axfather <- svyglm(father ~ sel, design = svydesign(ids = ~1, data = nlsy97_full[[1]], weights = ~weight))
#
#svyglm(paste0("afqt3 ~ hs +", paste(x, collapse = "+")), design = svydesign(ids = ~1, data = nlsy97_full[[1]], weights = ~weight))

# lm(afqt3 ~ hs +female+black+hispan+parinc+paredu+parast+livebio+father+afqt3+gpa+children18+substind+delinqind+rural97+south97+frcoll_75+frcoll_90+schstolen+schthreat+schfight, data = nlsy97_full[[1]] )
# 
# 0.0735572316
# 
# yxfather <- svyglm(comp ~ father + female + black + parinc+paredu + parast + rural97 + south97, design = svydesign(ids = ~1, data = nlsy97_full[[1]], weights = ~weight))
# 
# svyglm(father ~ sel +  female + black + parinc+paredu + parast + rural97 + south97, design = svydesign(ids = ~1, data = nlsy97_full[[1]], weights = ~weight))
# 
# # children by 18
# yxchildren18 <- svyglm(comp ~ children18, design = svydesign(ids = ~1, data = nlsy97_full[[1]], weights = ~weight))
# axchildren18 <- svyglm(children18 ~ sel, design = svydesign(ids = ~1, data = nlsy97_full[[1]], weights = ~weight))


###########################################
############ Draw contour plots ############
#############################################
my_size <- 3
u_name <- "Ability"

x <- seq(-1.1, 1.1 , by = 0.01) # consider difference in prevalence that's both
y <- seq(0, 2, by = 0.01) # consider effect of U on Y that's positive

# df <- expand.grid(x, y) %>%
#   `colnames<-`(c("x", "y")) %>%  
#   mutate(value= .14 + .57*x*y)
# 
# ggplot(df, aes(x, y)) +
#   geom_contour(aes(z = value), binwidth = .5,  colour = "black") +
#   metR::geom_text_contour(aes_(z = quote(value)), size = 3.5, colour = "black",  skip = 0, stroke = .5)

# tau0
df_t0 <- expand.grid(x, y) %>%
  `colnames<-`(c("x", "y")) %>%  
  mutate(value= overall_df[[8,2]] - x*y)

# tau1
df_t1 <- expand.grid(x, y) %>%
  `colnames<-`(c("x", "y")) %>%  
  mutate(value= overall_df[[9,2]] - x*y)

# tau2
df_t2 <- expand.grid(x, y) %>%
  `colnames<-`(c("x", "y")) %>%  
  mutate(value= overall_df[[10,2]] - x*y)

# tau3
df_t3 <- expand.grid(x, y) %>%
  `colnames<-`(c("x", "y")) %>%  
  mutate(value= overall_df[[11,2]] - x*y)

# delta0
df_d0 <- expand.grid(x, y) %>%
  `colnames<-`(c("x", "y")) %>%  
  mutate(value= overall_df[[1,2]] + overall_df[[5,2]]* x*y)

# delta1
df_d1 <- expand.grid(x, y) %>%
  `colnames<-`(c("x", "y")) %>%  
  mutate(value= overall_df[[2,2]] + overall_df[[6,2]]* x*y)

# delta2
df_d2 <- expand.grid(x, y) %>%
  `colnames<-`(c("x", "y")) %>%  
  mutate(value= overall_df[[3,2]] + overall_df[[7,2]]* x*y)

#### plots
row_number <- 21 # afqt3_bi

x1 <- as.numeric(out[[1]][row_number,3])
y1 <- as.numeric(out[[1]][row_number,2])
red <- sqrt(overall_df[[8,2]])

p_t0 <- ggplot(df_t0, aes(x, y)) +
  geom_contour(aes(z = value), binwidth = .5,  colour = "black") +
  metR::geom_text_contour(aes_(z = quote(value)), size = 3.5, colour = "black",  skip = 0, stroke = .5) +
  labs(x = "", y = "") +
  ggtitle(expression(paste(tau[0]))) +
  theme(axis.text.x = element_blank()) +
  annotate("point", x = x1, y = y1, color = "deepskyblue4", size = 2.5) +
  annotate("point", x = red, y = red, color = "deepskyblue4", size = 2.5) +
  annotate("text", x = red-0.05, y = red-0.12, label = paste0("(", round(red, digits = 2), ",", round(red, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") +
  annotate("text", x = x1 + 0.00, y = y1 + 0.25, label = paste0(u_name, " \n (", round(x1, digits = 2), ",", round(y1, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") 

x1 <- as.numeric(out[[2]][row_number,3])
y1 <- as.numeric(out[[2]][row_number,2])
red <- sqrt(overall_df[[9,2]])

p_t1 <- ggplot(df_t1, aes(x, y)) +
  geom_contour(aes(z = value), binwidth = .5,  colour = "black") +
  metR::geom_text_contour(aes_(z = quote(value)), size = 3.5, colour = "black",  skip = 0, stroke = .5) +
  labs(x = "", y = "") +
  ggtitle(expression(paste(tau[1]))) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  annotate("point", x = x1, y = y1, color = "deepskyblue4", size = 2.5) +
  annotate("point", x = red, y = red, color = "deepskyblue4", size = 2.5) +
  annotate("text", x = red-0.05, y = red-0.12, label = paste0("(", round(red, digits = 2), ",", round(red, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") +
  annotate("text", x = x1 + 0.01, y = y1 + 0.25, label = paste0(u_name, " \n(", round(x1, digits = 2), ",", round(y1, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") 


x1 <- as.numeric(out[[3]][row_number,3])
y1 <- as.numeric(out[[3]][row_number,2])
red <- sqrt(overall_df[[10,2]])

p_t2 <- ggplot(df_t2, aes(x, y)) +
  geom_contour(aes(z = value), binwidth = .5,  colour = "black") +
  metR::geom_text_contour(aes_(z = quote(value)), size = 3.5, colour = "black",  skip = 0, stroke = .5) +
  labs(x = "", y = "") +
  ylim(0,2)+
  ggtitle(expression(paste(tau[2]))) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  annotate("point", x = x1, y = y1, color = "deepskyblue4", size = 2.5) +
  annotate("point", x = red, y = red, color = "deepskyblue4", size = 2.5) +
  annotate("text", x = red-0.05, y = red-0.12, label = paste0("(", round(red, digits = 2), ",", round(red, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") +
  annotate("text", x = x1 + 0.00, y = y1 + 0.25, label = paste0(u_name, " \n(", round(x1, digits = 2), ",", round(y1, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") 

x1 <- as.numeric(out[[4]][row_number,3])
y1 <- as.numeric(out[[4]][row_number,2])
red <- sqrt(overall_df[[11,2]])

p_t3 <- ggplot(df_t3, aes(x, y)) +
  geom_contour(aes(z = value), binwidth = .5,  colour = "black") +
  metR::geom_text_contour(aes_(z = quote(value)), size = 3.5, colour = "black",  skip = 0, stroke = .5) +
  labs(x = "", y = "") +
  ggtitle(expression(paste(tau[3]))) +
  theme(axis.text.y = element_blank()) +
  annotate("point", x = x1, y = y1, color = "deepskyblue4", size = 2.5) +
  annotate("point", x = red, y = red, color = "deepskyblue4", size = 2.5) +
  annotate("text", x = red+0.198, y = red+0.035, label = paste0("(", round(red, digits = 2), ",", round(red, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") +
  annotate("text", x = x1 - 0.015, y = y1 + 0.25, label = paste0(u_name, " \n(", round(x1, digits = 2), ",", round(y1, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") 

## note for delta bias correction, no longer symmetric as one value has to be be negatie - i'll assume one is 
# overall_df[[1,2]] + overall_df[[5,2]]* x*y = 0 -> u_k = sqrt(-overall_df[[1,2]] +overall_df[[5,2]]) times -1 , factoring out negative 1


x1 <- as.numeric(out[[2]][row_number,3])
y1 <- as.numeric(out[[2]][row_number,2])
red <- sqrt(overall_df[[1,2]] / overall_df[[5,2]])

p_d0 <- ggplot(df_d0, aes(x, y)) +
  geom_contour(aes(z = value), binwidth = .2,  colour = "black") +
  metR::geom_text_contour(aes_(z = quote(value)), size = 3.5, colour = "black",  skip = 0, stroke = .5) +
  labs(x = "", y = "") +
  ylim(0,2)+
  ggtitle(expression(paste(Delta[0]))) +
  theme(axis.text.x = element_blank()) +
  annotate("point", x = x1, y = y1, color = "deepskyblue4", size = 2.5) +
  annotate("point", x = -red, y = red, color = "deepskyblue4", size = 2.5) +
  annotate("text", x = -red+0.12, y = red-0.12, label = paste0("(", round(red, digits = 2), ",", round(red, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") +
  annotate("text", x = x1 + 0.00, y = y1 + 0.25, label = paste0(u_name, " \n(", round(x1, digits = 2), ",", round(y1, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") 

x1 <- as.numeric(out[[3]][row_number,3])
y1 <- as.numeric(out[[3]][row_number,2])
red <- sqrt(overall_df[[2,2]] / overall_df[[6,2]])


p_d1 <- ggplot(df_d1, aes(x, y)) +
  geom_contour(aes(z = value), binwidth = .2,  colour = "black") +
  metR::geom_text_contour(aes_(z = quote(value)), size = 3.5, colour = "black",  skip = 0, stroke = .5) +
  labs(x = "", y = "") +
  ggtitle(expression(paste(Delta[1]))) +
  theme(axis.text.x = element_blank()) +
  annotate("point", x = x1, y = y1, color = "deepskyblue4", size = 2.5) +
  annotate("point", x = -red, y = red, color = "deepskyblue4", size = 2.5) +
  annotate("text", x = -red+0.12, y = red-0.12, label = paste0("(", round(red, digits = 2), ",", round(red, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") +
  annotate("text", x = x1 + 0.00, y = y1 + 0.25, label = paste0(u_name, " \n(", round(x1, digits = 2), ",", round(y1, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") 

x1 <- as.numeric(out[[4]][row_number,3])
y1 <- as.numeric(out[[4]][row_number,2])
red <- sqrt(overall_df[[3,2]] / overall_df[[7,2]])
# 1.367885 - oto extreme to fit on plot

p_d2 <- ggplot(df_d2, aes(x, y)) +
  geom_contour(aes(z = value), binwidth = .1,  colour = "black") +
  metR::geom_text_contour(aes_(z = quote(value)), size = 3.5, colour = "black",  skip = 0, stroke = .5) +
  labs(x = "", y = "") +
  ggtitle(expression(paste(Delta[2]))) +
  theme(axis.text.y = element_blank())+
  annotate("point", x = x1, y = y1, color = "deepskyblue4", size = 2.5) +
  # since 1.367885 - is too extreme to fit on plot
  # annotate("point", x = -red, y = red, color = "black", size = 2.5) +
  # annotate("text", x = -red+0.12, y = red-0.12, label = paste0("(", round(red, digits = 2), ",", round(red, digits = 2), ")"), size = my_size) +
  annotate("text", x = x1 + 0.00, y = y1 + 0.25, label = paste0(u_name, " \n(", round(x1, digits = 2), ",", round(y1, digits = 2), ")"), size = my_size,
           col = "deepskyblue4") 

left_caption <- expression(paste("Effect of ", U[k], " on ", Y, " (", beta[k], ")"))
bottom_caption <- expression(paste("Difference in Prevalance of ", U[k], " (", alpha[k], ")"))

cap0 = expression(paste("(", alpha[1], ", ", beta[1], " \n ) \n\n \n \n \n \n \n  f")) 

# Create a grid layout with the desired arrangement
sens_out = grid.arrange(
  arrangeGrob(p_t0, nullGrob(), widths = unit(c(1, 1), "null")),
  arrangeGrob(p_d0, p_t1, ncol = 2),
  arrangeGrob(p_d1, p_t2, ncol = 2),
  arrangeGrob(p_d2, p_t3, ncol = 2),
  ncol = 1
) %>%
  ggpubr::annotate_figure(bottom = textGrob(bottom_caption, gp = gpar(fontsize = 18)), 
                          left = textGrob(left_caption, rot = 90, gp = gpar(fontsize = 18)))
#right = textGrob("k=0\n\n\n\n\nk=1\n\n\n\n\nk=2\n\n\n\n\nk=3" , rot = 00, gp = gpar(fontsize = 10)))

## now removing tau 3

# sens_out = grid.arrange(
#   arrangeGrob(p_t0, p_d0, ncol = 2),
#   arrangeGrob(p_d1, p_t1, ncol = 2),
#   arrangeGrob(p_d2, p_t2, ncol = 2),
#   ncol = 1
# ) %>%
#   ggpubr::annotate_figure(bottom = textGrob(bottom_caption, gp = gpar(fontsize = 18)), 
#                           left = textGrob(left_caption, rot = 90, gp = gpar(fontsize = 18)))


sens_out

ggsave(sens_out, file = "Output/figsensitivity.png", width = 9.5, height = 11.8, bg = "white")




# ggplot(df, aes(x, y)) +
#   geom_contour(aes(z = value), binwidth = .025,  colour = "black") +
#   metR::geom_text_contour(aes_(z = quote(value)), size = 3.5, colour = "black",  skip = 0, stroke = .5) +
#   labs(x = xlabel, y = ylabel) +
#   geom_point(aes(x = x1, y = y1), col = "deepskyblue4") +
#   geom_point(aes(x = x2, y = y2), col = "deepskyblue4") +
#   geom_point(aes(x = 0.38, y = 0.38)) +
#   geom_text(aes(x = 0.52, y = 0.38, label = "(0.38,0.38)"), col = "deepskyblue4", size = 4) +
#   geom_text(aes(x = x1 + 0., y = y1 - 0.025, label = paste0("Urban Residence (", round(x1, digits = 2), ", ", round(y1, digits = 2), ")")),
#             col = "deepskyblue4", size = 4) +
#   geom_text(aes(x = x2, y = y2 + 0.025, label = paste0("Peers' College Expectation (", round(x2, digits = 2), ", ", round(y2, digits = 2), ")")),
#             col = "deepskyblue4", size = 4) +
#   ()

#ggsave(file = "figsensitivity2.png", width = 9, height = 6, bg = "white")




