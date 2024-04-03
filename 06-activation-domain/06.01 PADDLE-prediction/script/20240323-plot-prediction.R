# title: plot PADDLE prediction
# author: bin he
# date: 2024-03-23

# load libraries
require(tidyverse)
require(cowplot)

# read in data
datSc <- read_csv("../output/20240316-ScerPho4-PADDLE-prediction.csv")
datCg <- read_csv("../output/20240316-CglaPho4-PADDLE-prediction.csv")

# plot
p1 <- datSc %>% 
  ggplot(aes(x = Pos, y = Z_score)) +
  geom_hline(yintercept = c(4, 6), linetype = 2, linewidth = .8, color = "gray40") +
  geom_line(color = "steelblue", linewidth = 1.5) +
  scale_y_continuous(breaks = c(0, 4, 6, 10), limits = c(-0.25, 10)) +
  scale_x_continuous(breaks = c(0, 100, 200, 300), limits = c(0, 312), 
                     expand = expansion(mult = c(0, 0))) +
  xlab(NULL) + ylab("Z-score") +
  theme_cowplot() +
  theme(
    axis.text = element_text(face = 2),
    axis.title = element_text(face = 2)
  )
ggsave(filename = "../output/img/20240323-ScerPho4-PADDLE-for-R35.png", p1,
       width = 4, height = 1.7)
ggsave(filename = "../output/img/20240401-ScerPho4-PADDLE-for-paper.png", 
       p1 + panel_border(color = "gray10", size = 1) + 
         theme(axis.line = element_blank()),
       width = 7.3, height = 2.5)

p2 <- datCg %>% 
  ggplot(aes(x = Pos, y = Z_score)) +
  geom_hline(yintercept = c(4, 6), linetype = 2, linewidth = .8, color = "gray40") +
  geom_line(color = "steelblue", linewidth = 1.5) +
  scale_y_continuous(breaks = c(0, 4, 6, 10), limits = c(-0.25, 10)) +
  scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500), limits = c(0, 535),
                     expand = expansion(mult = c(0,0))) +
  xlab(NULL) + ylab("Z-score") +
  theme_cowplot() +
  theme(
    axis.text = element_text(face = 2),
    axis.title = element_text(face = 2)
  )
ggsave(filename = "../output/img/20240323-CglaPho4-PADDLE-for-R35.png", p2,
       width = 4, height = 1.7)
ggsave(filename = "../output/img/20240401-CglaPho4-PADDLE-for-paper.png",
       p2 + panel_border(color = "gray10", size = 1) + 
         theme(axis.line = element_blank()),
       width = 7.3, height = 2.5)

# compare multiple orthologs
my_read_paddle <- function(name){
  # read the file
  tmp <- read_csv(name, col_types = cols())
  # extend the scores to the two ends
  tmp <- tmp[c(rep(1, 26), 1:nrow(tmp), rep(nrow(tmp), 26)),] %>% 
    mutate(Pos = 1:n(), rel_pos = Pos/n())
}
files = dir("../output/", pattern = "20240316*")
names(files) <- str_split(files, pattern = "-", simplify = TRUE)[,2]
dat <- map_dfr(files, \(f) my_read_paddle(paste0("../output/", f)), .id = "Pho4")

# plot
bHLH <- tribble(
  ~Pho4, ~Feature, ~Begin, ~End,
  "ScerPho4", "bHLH", 250, 306,
  "SmikPho4", "bHLH", 249, 305,
  "CglaPho4", "bHLH", 470, 524,
  "CbraPho4", "bHLH", 424, 478,
  "LkluPho4", "bHLH", 419, 473,
  "CalbPho4", "bHLH", 594, 648
)
P2ID <- tribble(
  ~Pho4, ~Feature, ~Begin, ~End,
  "ScerPho4", "P2ID", 203, 249,
  "SmikPho4", "P2ID", 203, 248,
  "CglaPho4", "P2ID", 326, 469,
  "CbraPho4", "P2ID", 295, 423,
  "LkluPho4", "P2ID", 347, 418,
  "CalbPho4", "P2ID", 534, 593
) %>% mutate(Len = End - Begin + 1)
h = 0.5
act.colors <- c("gray80", "royalblue", "royalblue3")
names(act.colors) <- c("no act", "medium", "strong")
select.sps <- c("ScerPho4", "SmikPho4", "CglaPho4", "CbraPho4", "LkluPho4", "CalbPho4")
p3 <- dat %>% 
  filter(Pho4 %in% select.sps) %>% 
  mutate(Pho4 = factor(Pho4, levels = rev(select.sps)),
         paddle = cut(Z_score, breaks = c(-0.5, 4, 6, 12), 
                      labels = c("no act", "medium", "strong"))) %>% 
  ggplot(aes(x = Pos, y = Pho4, height = h)) +
  geom_tile(aes(fill = paddle), width = 2, linewidth = 0.5) +
  geom_tile(aes(x = (Begin + End)/2, y = Pho4, width = (End - Begin + 1)), 
            data = bHLH, height = h, fill = "yellowgreen") +
  scale_fill_manual("Activation", values = act.colors) +
  theme_cowplot() +
  theme(
    axis.line = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )
ggsave("../output/img/20240324-diverse-Pho4-PADDLE-for-R35.png", p3, width = 4, height = 3)  
