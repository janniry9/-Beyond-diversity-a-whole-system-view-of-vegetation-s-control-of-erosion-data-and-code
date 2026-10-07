#Beyond diversity: a whole system view of vegetation’s control of erosion

#Belen, J.C; Livsey, D.N; Martin, C.B; Groh, T.A; Cardinale, B.J

####BOOTSTRAPPING R code for computing the RY, from "resampling" of plant individuals. All variables created are directly taken from the Loreau 2001 paper ####

#data set-up#
#load in data set#
relative_yield_1 <- read_csv("data/relative_yield_1.csv")
View(relative_yield_1)

n_unit <- 18   # plants per unit #
B <- 1000 #bootstrap iterations #

#monoculture data#
mono_indiv <- relative_yield %>%
  filter(richness == 1) %>%
  mutate(root_biomass = as.numeric(root_biomass)) %>%
  filter(is.finite(root_biomass), root_biomass > 0) %>%
  group_by(treatment) %>%
  summarise(
    indiv_vals = list(root_biomass),
    M_point = mean(root_biomass) * n_unit,   
    .groups = "drop"
  )


# polyculture data: creation of O and p#
polycultures <- relative_yield %>%
  filter(richness > 1) %>%
  group_by(unit_id) %>%
  mutate(N_total = n()) %>%
  group_by(unit_id, richness, treatment) %>%
  summarise(
    O = sum(root_biomass, na.rm = TRUE),
    Ni = n(),
    N_total = first(N_total),
    p = Ni / N_total,
    .groups = "drop"
  )


#Bootstrapping function: creation of E, RY, D #

boot_species_metrics_indiv <- function(O, p, indiv_vals, n_unit = 18, B = 1000) {
  vals <- as.numeric(unlist(indiv_vals))
  vals <- vals[is.finite(vals) & vals > 0]
  
  
  M_draw <- replicate(B, sum(sample(vals, size = n_unit, replace = TRUE))) 
  
  E_draw  <- p * M_draw
  RY_draw <- O / M_draw
  D_draw  <- (O - E_draw) / E_draw
  
  ok <- is.finite(M_draw) & M_draw > 0 &
    is.finite(E_draw) & E_draw != 0 &
    is.finite(RY_draw) &
    is.finite(D_draw)
  
  tibble(
    RY_mean = mean(RY_draw[ok], na.rm = TRUE),
    RY_sd   = sd(RY_draw[ok],   na.rm = TRUE),
    D_mean  = mean(D_draw[ok],  na.rm = TRUE),
    D_sd    = sd(D_draw[ok],    na.rm = TRUE),
    E_mean  = mean(E_draw[ok],  na.rm = TRUE),
    E_sd    = sd(E_draw[ok],    na.rm = TRUE)
  )
}

# creating: set points of E, RY, D #

set.seed(123)

species_table <- polycultures %>%
  left_join(mono_indiv, by = "treatment") %>%
  filter(!map_lgl(indiv_vals, ~ is.null(.x) || length(.x) == 0)) %>%
  mutate(
    # point estimates
    M_point = M_point,
    E  = p * M_point,
    RY = O / M_point,
    D  = (O - E) / E,
    
    # bootstrap summaries using individual pools
    stats = pmap(list(O = O, p = p, indiv_vals = indiv_vals),
                 ~ boot_species_metrics_indiv(..1, ..2, ..3, n_unit = n_unit, B = B))
  ) %>%
  unnest(stats)


#unit level variables (RYT, OT, ET, DT) #

boot_unit_metrics_indiv <- function(df_unit, B = 1000, n_unit = 18) {
  Oi <- as.numeric(df_unit$O)
  pi <- as.numeric(df_unit$p)
  pools <- df_unit$indiv_vals
  
  OT_point  <- sum(Oi, na.rm = TRUE)
  ET_point  <- sum(df_unit$E,  na.rm = TRUE)
  RYT_point <- sum(df_unit$RY, na.rm = TRUE)
  DT_point  <- (OT_point - ET_point) / ET_point
  

  M_mat <- do.call(cbind, lapply(pools, function(v) {
    vals <- as.numeric(unlist(v))
    vals <- vals[is.finite(vals) & vals > 0]
    replicate(B, sum(sample(vals, size = n_unit, replace = TRUE)))
  }))
  M_mat <- as.matrix(M_mat)
  storage.mode(M_mat) <- "double"
  
  ET_draw  <- rowSums(sweep(M_mat, 2, pi, `*`), na.rm = TRUE)
  RYT_draw <- rowSums(sweep(1 / M_mat, 2, Oi, `*`), na.rm = TRUE)
  DT_draw  <- (OT_point - ET_draw) / ET_draw
  
  ok <- is.finite(ET_draw) & ET_draw != 0 &
    is.finite(RYT_draw) &
    is.finite(DT_draw)
  
  tibble(
    OT  = OT_point,
    ET  = ET_point,
    RYT = RYT_point,
    DT  = DT_point,
    
    ET_mean  = mean(ET_draw[ok],  na.rm = TRUE),
    ET_sd    = sd(ET_draw[ok],    na.rm = TRUE),
    RYT_mean = mean(RYT_draw[ok], na.rm = TRUE),
    RYT_sd   = sd(RYT_draw[ok],   na.rm = TRUE),
    DT_mean  = mean(DT_draw[ok],  na.rm = TRUE),
    DT_sd    = sd(DT_draw[ok],    na.rm = TRUE)
  )
}

unit_table_boot <- species_table %>%
  group_by(unit_id, richness) %>%
  group_modify(~ boot_unit_metrics_indiv(.x, B = B, n_unit = n_unit)) %>%
  ungroup()

# join the two tables and write the csv #

final_table <- species_table %>%
  left_join(unit_table_boot, by = c("unit_id", "richness")) %>%
  arrange(unit_id, treatment)

final_table_csv <- final_table %>%
  select(-any_of(c("indiv_vals"))) %>%
  mutate(across(where(is.list), ~ NULL))

write.csv(final_table_csv, "RY_bootstrap_indiv.csv", row.names = FALSE)

