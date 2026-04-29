library(mgcv)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)
library(data.table)

apply_capping <- function(df, ref_df) {
  num_cols <- names(ref_df)[sapply(ref_df, is.numeric)]
  num_cols <- num_cols[!num_cols %in% c("Total_Cases", "year")]
  for (col in num_cols) {
    if (col %in% names(df)) {
      min_val <- min(ref_df[[col]], na.rm = TRUE)
      max_val <- max(ref_df[[col]], na.rm = TRUE)
      df[[col]] <- pmax(min_val, pmin(max_val, df[[col]]))
    }
  }
  return(df)
}

get_predictions_with_ci <- function(model, newdata, scenario_name) {
  preds <- predict(model, newdata, type = "link", se.fit = TRUE)
  newdata$link_fit <- preds$fit
  newdata$link_se <- preds$se.fit
  newdata$cases <- exp(newdata$link_fit)
  newdata$lower <- exp(newdata$link_fit - 1.96 * newdata$link_se)
  newdata$upper <- exp(newdata$link_fit + 1.96 * newdata$link_se)
  newdata$Scenario <- scenario_name
  return(newdata)
}

cholera_data$country <- as.factor(cholera_data$country)
cholera_data$Total_Cases <- as.integer(cholera_data$Total_Cases)

refined_model <- gam(
  Total_Cases ~ s(GDP_per_capita, k = 5) + 
    s(basic_drinking_water_access_pct, k = 5) + 
    s(open_defecation_percent, k = 5) + 
    s(precipitation, k = 5) + 
    s(extreme_temp_days, k = 5) +
    ti(precipitation, open_defecation_percent) + 
    s(country, bs = "re"),
  family = nb(),
  data = cholera_data,
  method = "REML"
)

years_future <- 2025:2054
scenarios <- c("BAU", "HTC", "INT")
projection_results <- data.frame()

for (yr in years_future) {
  base_step <- future_base_data %>% filter(year == yr)
  
  data_bau <- base_step %>%
    mutate(Scenario = "BAU") %>%
    apply_capping(cholera_data)
  
  res_bau <- get_predictions_with_ci(refined_model, data_bau, "BAU")
  projection_results <- rbind(projection_results, res_bau)
  
  data_htc <- base_step %>%
    mutate(
      basic_drinking_water_access_pct = basic_drinking_water_access_pct + (yr - 2024) * 0.5,
      open_defecation_percent = open_defecation_percent - (yr - 2024) * 0.5,
      Scenario = "HTC"
    ) %>%
    mutate(
      basic_drinking_water_access_pct = pmin(100, basic_drinking_water_access_pct),
      open_defecation_percent = pmax(0, open_defecation_percent)
    ) %>%
    apply_capping(cholera_data)
  
  res_htc <- get_predictions_with_ci(refined_model, data_htc, "HTC")
  projection_results <- rbind(projection_results, res_htc)
  
  data_int <- base_step %>%
    mutate(
      basic_drinking_water_access_pct = basic_drinking_water_access_pct + (yr - 2024) * 1.0,
      open_defecation_percent = open_defecation_percent - (yr - 2024) * 1.0,
      Scenario = "INT"
    ) %>%
    mutate(
      basic_drinking_water_access_pct = pmin(100, basic_drinking_water_access_pct),
      open_defecation_percent = pmax(0, open_defecation_percent)
    ) %>%
    apply_capping(cholera_data)
  
  res_int <- get_predictions_with_ci(refined_model, data_int, "INT")
  projection_results <- rbind(projection_results, res_int)
}

annual_summary <- projection_results %>%
  group_by(year, Scenario) %>%
  summarise(
    total_cases = sum(cases, na.rm = TRUE),
    total_lower = sum(lower, na.rm = TRUE),
    total_upper = sum(upper, na.rm = TRUE),
    .groups = "drop"
  )

plot_main <- ggplot(annual_summary, aes(x = year, y = total_cases, color = Scenario, fill = Scenario)) +
  geom_ribbon(aes(ymin = total_lower, ymax = total_upper), alpha = 0.15, color = NA) +
  geom_line(size = 1.2) +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("BAU" = "#d95f02", "HTC" = "#7570b3", "INT" = "#1b9e77")) +
  scale_fill_manual(values = c("BAU" = "#d95f02", "HTC" = "#7570b3", "INT" = "#1b9e77")) +
  labs(title = "Long-term Cholera Projections (2025-2054)", x = "Year", y = "Annual Total Cases") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

baseline_2024_val <- sum(cholera_data$Total_Cases[cholera_data$year == 2024], na.rm = TRUE)
target_elimination <- baseline_2024_val * 0.1

rates_test <- seq(0, 15, by = 0.5)
target_seeking_data <- data.frame()

data_2030_raw <- future_base_data %>% filter(year == 2030)

for (r in rates_test) {
  sim_2030 <- data_2030_raw %>%
    mutate(
      basic_drinking_water_access_pct = pmin(100, basic_drinking_water_access_pct + (6 * r)),
      open_defecation_percent = pmax(0, open_defecation_percent - (6 * r))
    ) %>%
    apply_capping(cholera_data)
  
  pred_2030 <- get_predictions_with_ci(refined_model, sim_2030, paste0(r, "%"))
  
  target_seeking_data <- rbind(target_seeking_data, data.frame(
    Acceleration_Rate = r,
    Total_Cases = sum(pred_2030$cases),
    Lower = sum(pred_2030$lower),
    Upper = sum(pred_2030$upper)
  ))
}

required_rate <- target_seeking_data$Acceleration_Rate[which(target_seeking_data$Total_Cases <= target_elimination)[1]]

plot_target <- ggplot(target_seeking_data, aes(x = Acceleration_Rate, y = Total_Cases)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#1b9e77", alpha = 0.2) +
  geom_line(color = "#1b9e77", size = 1.5) +
  geom_hline(yintercept = target_elimination, linetype = "dashed", color = "#e7298a", size = 1) +
  geom_vline(xintercept = required_rate, linetype = "dotted", color = "black") +
  annotate("text", x = 12, y = target_elimination * 1.5, label = "90% Reduction Target", color = "#e7298a") +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Target-Seeking Analysis for 2030 Goals",
    x = "Annual WASH Acceleration Rate (%)",
    y = "Projected Cases in 2030"
  ) +
  theme_minimal(base_size = 14)

regional_analysis <- projection_results %>%
  filter(year %in% c(2030, 2040, 2050)) %>%
  group_by(country, year, Scenario) %>%
  summarise(cases = sum(cases), .groups = "drop")

final_table <- annual_summary %>%
  mutate(
    Period = case_when(
      year <= 2034 ~ "2025-2034",
      year <= 2044 ~ "2035-2044",
      TRUE ~ "2045-2054"
    )
  ) %>%
  group_by(Period, Scenario) %>%
  summarise(
    Mean_Cases = mean(total_cases),
    Min_Cases = min(total_cases),
    Max_Cases = max(total_cases),
    .groups = "drop"
  )

print(plot_main)
print(plot_target)
write.csv(annual_summary, "annual_projections.csv", row.names = FALSE)
write.csv(target_seeking_data, "target_seeking_results.csv", row.names = FALSE)
write.csv(final_table, "decadal_summary_table.csv", row.names = FALSE)

cat("Analysis Complete. Required acceleration rate to hit 2030 target:", required_rate, "%\n")
