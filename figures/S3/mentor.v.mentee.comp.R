library(coda)
library(ggplot2)
library(ggrepel)

files.student <- list.files("../results/exponential.prior/trainee_results", 
                            full.names = TRUE)
files.mentor <- list.files("../results/exponential.prior/mentor_results", 
                           full.names = TRUE)

# Initialize storage - ADDED mentor.prob and student.prob columns
results <- data.frame(
  clade = character(length(files.student)),
  mentor.lower = numeric(length(files.student)),
  mentor.upper = numeric(length(files.student)),
  mentor.median = numeric(length(files.student)),
  mentor.prob = numeric(length(files.student)),   # New column
  student.lower = numeric(length(files.student)),
  student.upper = numeric(length(files.student)),
  student.median = numeric(length(files.student)),
  student.prob = numeric(length(files.student))   # New column
)

for (i in 1:length(files.student)) {
  
  # Get clade name
  x <- gsub("\\.rds$", "", basename(files.student[i]), ignore.case = TRUE)
  results$clade[i] <- substr(x, 1, nchar(x) - 9)
  # Process mentor file
  res.mentor <- read.csv(files.mentor[i])[-1:-100,]
  total.mentor <- res.mentor$asc1 + res.mentor$desc1
  hpd.mentor <- HPDinterval(as.mcmc(total.mentor), prob = 0.95)
  results$mentor.lower[i] <- hpd.mentor[1, "lower"]
  results$mentor.upper[i] <- hpd.mentor[1, "upper"]
  results$mentor.median[i] <- median(total.mentor)
  results$mentor.prob[i] <- mean(res.mentor$p)      # Calculate mean probability
  
  # Process student file
  res.student <- read.csv(files.student[i])[-1:-100,]
  total.student <- res.student$asc1 + res.student$desc1
  hpd.student <- HPDinterval(as.mcmc(total.student), prob = 0.95)
  results$student.lower[i] <- hpd.student[1, "lower"]
  results$student.upper[i] <- hpd.student[1, "upper"]
  results$student.median[i] <- median(total.student)
  results$student.prob[i] <- mean(res.student$p)    # Calculate mean probability
}

# LOGIC:
# 1. Identify statistical differences (no overlap)
results$is_different <- with(results, student.lower > mentor.upper | student.upper < mentor.lower)

# 2. Highlight if different OR if it is "Solanaceae"
results$highlight <- results$is_different | grepl("Solanaceae", results$clade, ignore.case = TRUE)

# Plot
ggplot(results, aes(x = mentor.median, y = student.median)) +
  # 1. Diagonal reference line (Bottom Layer)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  
  # 2. Points (Middle Layer - Now plotted BEFORE lines so lines sit on top)
  geom_point(aes(alpha = highlight, 
                 color = I(ifelse(highlight, "black", "gray40")))) +
  
  # 3. MENTOR Error Bars (Horizontal) - Plotted ON TOP
  geom_segment(aes(x = mentor.lower, xend = mentor.upper, 
                   y = student.median, yend = student.median,
                   alpha = highlight,
                   color = I(ifelse(highlight, "red", "gray40"))),
               lwd = 0.8) +
  
  # 4. STUDENT Error Bars (Vertical) - Plotted ON TOP
  geom_segment(aes(x = mentor.median, xend = mentor.median, 
                   y = student.lower, yend = student.upper,
                   alpha = highlight,
                   color = I(ifelse(highlight, "darkblue", "gray40"))),
               lwd = 0.8) +
  
  # Alpha Scaling: Non-highlighted = 0.5, Highlighted = 1.0
  scale_alpha_manual(values = c(`FALSE` = 0.5, `TRUE` = 1), guide = "none") +
  
  labs(x = "Mentor: Dysploidy Rates (95% HPD; units of tree length)",
       y = "Student: Sum of Changes (95% HPD; units of tree length)") +
  theme_bw()
