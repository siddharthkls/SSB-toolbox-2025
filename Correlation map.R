# Load required libraries
library(ggplot2)
library(reshape2)

# Read the CSV file
data <- read.csv("D:/Toolbox/Result2.csv", header = TRUE)

# Rename columns for consistency
colnames(data) <- c("Sector", "Mean_Conservation", "Mean_Correlation", "Negative_Correlation", "Positive_Correlation",
                    "Corrected_sectors", "Sectors_1", "Negative_correlation_1","Positive_correlation_1", "Corrected_sectors1")

# Select relevant columns
cor_data <- data[, c("Corrected_sectors", "Negative_Correlation","Positive_Correlation")]

# Melt to long format
long_data <- melt(cor_data, id.vars = "Corrected_sectors", variable.name = "Type", value.name = "Percentage")

# Remove NA values
long_data <- na.omit(long_data)

# Rename factor levels in 'Type' for clarity
long_data$Type <- factor(long_data$Type,
                         levels = c("Positive_Correlation", "Negative_Correlation"),
                         labels = c("Positive", "Negative"))

# Rename Sector column
colnames(long_data)[colnames(long_data) == "Corrected_sectors"] <- "Sector"
long_data$Sector <- as.factor(long_data$Sector)

# Save to PNG
png("D:/Toolbox/Result2_1500.png", width = 1600, height = 1000, res = 150)

# Plot
ggplot(long_data, aes(x = Sector, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", color = "gray60", linewidth = 0.4) +
  scale_fill_manual(values = c("Positive" = "lightblue", "Negative" = "brown")) +
  labs(title = "Percentage of Positive and Negative Correlations per Sector",
       x = "Sector",
       y = "Percentage (%)",
       fill = "Correlation Type") +
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.background = element_blank(),
    legend.key = element_blank()
  )

dev.off()
