setwd("/Users/bfj994/Documents/barcodeMiner/")

library(ggplot2)

# Read the .prof file
prof <- read.delim("dhigh5.dat", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "")

# Extract C>T values (removing "[0..0]" part)
ct_freq <- as.numeric(sub(" .*", "", prof$`C>T`))
positions <- prof[[1]]

# Extend to position 20 by repeating the 6th value
last_value <- ct_freq[7]
extended_positions <- c(positions, 7:20)
extended_ct <- c(ct_freq, rep(last_value, length(7:20)))

# Create a data frame for ggplot
df <- data.frame(
  Position = extended_positions,
  Frequency = extended_ct
)

# Plot using ggplot2
ggplot(df, aes(x = Position, y = Frequency)) +
  geom_line(color = "#D55E00", size = 1) +
  geom_point(color = "#D55E00", size = 2) +
  labs(
    title = "5′ Damage (C->T)",
    x = "Position from 5′ end",
    y = "C>T Frequency"
  ) +
  theme_bw(base_size = 14)

# Read the 3′ end damage profile (allow flexible whitespace separation)
prof3 <- read.delim("dhigh3.dat", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "")

# Extract G>A column (adjust if your column is named differently)
ga_freq <- as.numeric(sub(" .*", "", prof3$`G>A`))

# Extract the position column (assumed to be the first column)
positions <- prof3[[1]]

# Find index of position -6
pos_index <- which(positions == -6)
ga_value_at_neg6 <- ga_freq[pos_index]

# Create extended positions from current minimum to -20
extended_positions <- seq(min(positions), -20, by = -1)

# Fill values: keep existing, then extend using value at -6
missing_count <- length(extended_positions) - length(positions)
extended_ga <- c(ga_freq, rep(ga_value_at_neg6, missing_count))

# Create a data frame for ggplot
df3 <- data.frame(
  Position = extended_positions,
  Frequency = extended_ga
)

# Plot using ggplot2
ggplot(df3, aes(x = Position, y = Frequency)) +
  geom_line(color = "#0072B2", size = 1) +
  geom_point(color = "#0072B2", size = 2) +
  labs(
    title = "3′ End DNA Damage Curve (G>A)",
    x = "Position from 3′ end",
    y = "G>A Frequency"
  ) +
  theme_minimal(base_size = 14)


########## Frag length distirbution 
# Read your fragment lengths
frag_lengths <- read.table("fragLen.txt", header = F)
colnames(frag_lengths) <- c("Count", "Length")
str(frag_lengths)
# Plot the fragment length distribution
ggplot(frag_lengths, aes(x = Length, y = Count)) +
  geom_col(fill = "#0072B2") +
  labs(
    title = "Fragment Length Distribution",
    x = "Fragment Length (bp)",
    y = "Count"
  ) +
  theme_bw(base_size = 14)
