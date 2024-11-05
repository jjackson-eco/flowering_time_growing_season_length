#####################################################
##                                                 ##
##            Optimal Flowering Time               ##
##                                                 ##
##               Theory Schematics                 ##
##                                                 ##
##              JP & JJ - July 2024                ##
##                                                 ##
#####################################################
rm(list = ls())
options(width = 100)

library(ggplot2)
library(egg)
library(ggforce)
library(ggpattern)
library(latex2exp) 

#_______________________________________________________________________________
#### 1. f variation ####

# GSL gradient
T_seq <- seq(70, 160, length.out = 100)

# shape parameters for f(T); just varying "b"
a <- 0.02
b_1 <- 0 # this equates to constant f i.e linear
b_2 <- -0.04
c <- -160 # "anchor point", 

# three variants of the f(T) function
f_1 <- a
f_2 <- a/(1 + exp(-b_2*(T_seq + c))) # curved f(T)

# three variants of Flowering Time function
t2_1 <- T_seq - 1/(f_1)
t2_2 <- T_seq - 1/(f_2)

inset_f_plot <- ggplot()+
  geom_line(aes(x = T_seq, y = f_1), size = 1, linetype = "dashed", col = "gray50") +
  geom_line(aes(x = T_seq, y = f_2), size = 1, col = "black") +
  labs(x = expression(italic(T)), y = expression(italic("f"))) +
  ylim(0.003, a + 0.007) +
  annotate('text', label = TeX('$f=a$'), x = 75, y = 0.025, size = 4.5, hjust = 0, col = "gray50") +
  annotate('text', label = TeX(r"($f(T) = \frac{a}{1 + e^{-b(T + c)}})"), parse = T, x = 75, y = 0.009, size = 4, hjust = 0, col = "black") +
  theme_bw() +
  theme(plot.background = element_rect(fill = "grey90", colour = "black"),
        panel.background = element_rect(fill = "grey90"), 
        panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(colour = "black", size = 12)
  )

variable_f_t2 <- ggplot()+
  geom_line(aes(x = T_seq, y = t2_1), size = 1.5, linetype = "dashed", col = "gray50") +
  geom_line(aes(x = T_seq, y = t2_2), size = 1.5, col = "black") +
  annotation_custom(
    ggplotGrob(inset_f_plot), 
    xmin = 67, xmax = 108, ymin = 63, ymax = 109.5) +
  labs(x = expression(paste("Growing Season Length (", italic("T"), ") [days]")),
       y = expression(paste("Flowering time from start of season (", italic(t[italic("2")]), ")"))) +
  annotate('text', label = "t[2]==T-1/f", parse = T, x = 130, y = 108, size = 6, hjust = 0, col = "gray50") +
  annotate('text', label = "t[2]==T-1/f(T)", parse = T, x = 130, y = 55, size = 6, hjust = 0, col = "black") +
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        axis.text = element_text(colour = "black", size = 12),
        axis.title = element_text(colour = "black", size = 14))

ggsave(variable_f_t2, 
       filename = "output/variable_f_fig20240930.png",
       width = 15, height = 12, units = "cm", dpi = 700)

#_______________________________________________________________________________
#### 2. Iwasa-Cohen perennial schedule ####

T <- 200 #GSL

S0 <- 2 #storage size at beginning of season
S_T <- 20 #storage size at end of season
t1 <- 20 #first switching time; when storage use ends
t2 <- 110 #second switching time; "flowering time"

S_0_t1 <- S0 * - exp(seq(0, t1, by = 1) / t1) - S0*-exp(1)
S_t1_t2 <- rep(0, t2-t1)
S_t2_T <- seq(0, S_T, length.out = T-t2)
S_function <- c(S_0_t1, S_t1_t2, S_t2_T)

F_0_t1 <- max(S_0_t1) - S_0_t1
F_t1_t2 <- (c(seq(0, t2-t1-1, by = 1))^2) / 1000 + max(F_0_t1)
F_t2_T <- rep(max(F_t1_t2), T-t2)
F_function <- c(F_0_t1, F_t1_t2, F_t2_T)

schedule_schematic <- ggplot() +
  # S and F functions
  geom_line(aes(x = seq(0, T, by =1), y = S_function), col = "black", size = 1) +
  geom_line(aes(x = seq(0, T, by =1), y = F_function), col = "black", size = 1) +
  geom_ribbon(aes(x = seq(0, T, by = 1), ymin = 0, ymax = S_function), fill = "#71350B", alpha = 0.3) +
  geom_ribbon(aes(x = seq(0, T, by = 1), ymin = 0, ymax = F_function), fill = "#0D8922", alpha = 0.3) +
  
  # Labelling functions 
  annotate('text', label = expression(italic(S[n])), x = 167, y = 15, size = 5, col = "black") +
  annotate('text', label = expression(italic(F[n])), x = 65, y = 7, size = 5, col = "black") +
  
  # Years label
  annotate('text', label = "Year n",  x = 3, y = S_T, hjust = 0, fontface = "italic", size = 5, col = "grey30") +
  annotate('text', label = "Year n+1",  x = T + 3, y = S_T, hjust = 0, fontface = "italic", size = 5, col = "grey30") +
  
  # Year n+1 gray-out
  annotate("rect", xmin = T, xmax = Inf, ymin = 0, ymax = S_T + 1, fill = "grey90", alpha = 0.4) +  
  
  # S_n(0) point
  geom_point(aes(x = 0, y = max(S_0_t1)), size = 3) +
  annotate('text', label = expression(italic(S)[italic(n)]*italic("(0)")), x = 11, y = 4.5, size = 5, parse = TRUE) +

  # Highlighting switching points
  #annotate('segment', x = t1, xend = t1, y = 0, yend = S_T - 3, linetype = "dotted", size = 1) +
  annotate('segment', x = t2, xend = t2, y = 0, yend = S_T - 4, linetype = "dashed", size = 0.7, col = "black") +
  annotate('text', label = "Flowering Time", x = t2, y = S_T - 3, hjust = 0.5, size = 5, fontface = "bold", col = "black") +
  
  # T vertical boundary
  annotate('segment', x = T, xend = T, y = 0, yend = S_T + 1, size = 0.4) +
  
  # R and gamma 
  annotate('segment', x = T, xend = T, y = S_T, yend = S_T - 6, arrow = arrow(type = "closed", length = unit(2.5, "mm")), size = 1.25) + 
  annotate('text', label = expression(italic(R[n])), x = T + 5.5, y = S_T - 3, vjust = 0.5, size = 5) +
  annotate('segment', x = T, xend = T + 10, y = S_T - 6, yend = S_T - 12, arrow = arrow(type = "closed", length = unit(2.5, "mm")), size = 1.25) + 
  annotate('text', label = expression(italic(gamma[n])), x = T + 9.5, y = S_T - 8.5, vjust = 0.5, size = 5) +
  
  # S_[n+1](0) point
  geom_point(aes(x = T + 10.5, y = S_T - 12.5), size = 3) +
  annotate('text', label = expression(italic(S)[italic("n+1")]*italic("(0)")), x = T + 23, y = S_T - 11.5, size = 5, parse = T) +
  
  # S_[n+1] curve
  geom_curve(aes(x = T + 10, xend = T + 30, y = S_T - 12.5, yend = 3.5), curvature = -0.2, size = 0.7, linetype = "dotted") +
  
  # Axes and labels
  labs(title = "Optimal annual schedule of a perennial life history", 
       x = expression(paste("Time within growing season (", italic("t"), ")")),
       y = "Size") +
  scale_y_continuous(limits = c(-1.3, S_T + 1), expand = c(0, -1.3, 0, 0)) +
  scale_x_continuous(breaks = c(0, t1, t2, T, T+10), labels = c("0", expression(italic(t[1])), expression(italic(t[2])), expression(italic(T)), "0"), limits = c(0,T+35), expand = c(0, 0, 0, 0)) +
  #geom_point(aes(x = t1, y = -1.3), size = 9, pch = 1, col = "black") +
  geom_point(aes(x = t2, y = -1.3), size = 9, pch = 1, col = "black") +
  
  # Discontinuous x-axis symbol (hack method)
  annotate('text', label = "//", x = T + 5, y = 0, size = 4) +
  
  coord_cartesian(clip = "off") +
  
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 17, face = "bold"),
        axis.text.x = element_text(color = c("black", "black", "black", "black")),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.tag.position = c(0.96,0.03))

ggsave(schedule_schematic, 
       filename = "output/schedule_schematic_fig20240930.png",
       width = 20, height = 9, units = "cm", dpi = 700)
