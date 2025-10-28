install.packages("minpack.lm")
install.packages("ggplot2")
install.packages("zoo")
library(minpack.lm)
library(ggplot2)
library(zoo)

#selecting and formatting data
print("Choose the data you'd like to analyze:")
file_path <- file.choose()
abs_spectra <- read.csv(file_path, header=TRUE)
print(file_path)

wavelengths1=abs_spectra[, 1]
abs_vectors <- lapply(abs_spectra[, -1], as.vector)

#use this to determine which spectrometer was used to collect data
target_index=1
min_wavelength <- wavelengths1[target_index]
print(min_wavelength)

#loading reference library
print("Choose the corresponding reference spectra:")
ref_file_path <- file.choose()
references <- read.csv(ref_file_path, header=TRUE)
print(ref_file_path)

wavelengths2 <- references[, 1]
reference_vectors <- lapply(references[, -1], as.vector)

min_wavelength_ref <- wavelengths2[target_index]
print(min_wavelength_ref)

print("Do the starting wavelengths match?")
if (min_wavelength == min_wavelength_ref) {
  print("yes, proceed")
} else {
  warning("no, mismatch")
}

#plot cytochrome bH to confirm that the references are being loaded
cyt_bH=reference_vectors[[6]]
plot(wavelengths2, cyt_bH, type = "l", xlab = "Wavelength (nm)", ylab = "OD", main = "Cytochrome bH")

#can also use ggplot to plot multiple references
ggplot(references, aes(x = wavelengths2)) +
  geom_line(aes(y = reference_vectors[[6]], color = "Cytochrome bH")) +
  geom_line(aes(y = reference_vectors[[7]], color = "Cytochrome bL")) +
  labs(x = "Wavelength (nm)", y = "OD", color = "Reference") +
  scale_color_manual(values = c("Cytochrome bH" = "blue", "Cytochrome bL" = "red")) +
  theme_minimal()

#establish the bandwidth used for analyses
start_wavelength <- 540
end_wavelength <- 580

#plot absorbance over time
num_vectors <- length(abs_vectors)

#select the wavelength you want to evaluate over time
target_wavelength <- 550
target_index <- which.min(abs(wavelengths2 - target_wavelength))
print(target_index)

#extract intensities at the target wavelength from all spectra
intensities <- sapply(abs_vectors, function(vector) {
  #return the intensity at the target index
  vector[target_index]
})

# Create a sequence of time points (assuming 1 sec per spectrum)
time <- seq(0, num_vectors-1)

#plot intensities as a function of time
plot(time, intensities, type = "l", xlab = "Time (seconds)", ylab = "Absorbance", main = paste("Intensity at", target_wavelength, "nm over Time"))

#find the indices corresponding to the start and end wavelengths
start_index <- which.min(abs(wavelengths2 - start_wavelength))
end_index <- which.min(abs(wavelengths2 - end_wavelength))
print(start_index)
print(end_index)


#crop the reference spectra to the desired bandwidth
subset_reference_vectors <- list()

for (i in seq_along(reference_vectors)) {
  # Subset the current vector to include only the desired range
  subset_vector <- reference_vectors[[i]][start_index:end_index]
  # Add the subset to the list
  subset_reference_vectors[[i]] <- subset_vector
}

wavelengths3=wavelengths2[start_index:end_index]

plot(wavelengths3, subset_reference_vectors[[12]], type = "l", xlab = "Wavelength", ylab = "Absorbance", main = "Cytochrome c")

#validate peak wavelength of reference spectra
max_index_c <- which.max(subset_reference_vectors[[12]])
peak_wav_c=wavelengths3[max_index_c]
print(peak_wav_c)

max_index_bH=which.max(subset_reference_vectors[[6]])
peak_wav_bH=wavelengths3[max_index_bH]
print(peak_wav_bH)

max_index_bL=which.max(subset_reference_vectors[[7]])
peak_wav_bL=wavelengths3[max_index_bL]
print(peak_wav_bL)

#perform a rolling average before fitting
abs_matrix <- as.matrix(abs_spectra[, -1])

print("How many consecutive spectra do you want to average over time?")
window_width <- 5

#calculate the total number of spectra
num_spectra <- ncol(abs_matrix)
print(num_spectra)

#initialize a matrix to store the rolling average
rolling_avg <- matrix(NA, nrow = nrow(abs_matrix), ncol = num_spectra - (window_width - 1))

#take the rolling average using the specified width over all spectra
for (i in 1:(num_spectra - (window_width - 1))) {
  rolling_avg[, i] <- rowMeans(abs_matrix[, (i):(i + window_width - 1)], na.rm = TRUE)
}

#creating a list of vectors from the rolling average
exp_spec_vectors <- data.frame(rolling_avg)

print(paste("Wavelength Range:", start_wavelength, "-", end_wavelength, sep = ""))

#create empty list to store the new vectors
subset_exp_vectors <- list()

#creating a subset (by wavelength) of the AVERAGED experimental spectra over time
for (i in seq_along(exp_spec_vectors)) {
  subset_exp_vector <- exp_spec_vectors[[i]][start_index:end_index]
  subset_exp_vectors[[i]] <- subset_exp_vector
}

plot(wavelengths3, subset_exp_vectors[[12]], type = "l", xlab = "Wavelength", ylab = "Intensity", main = "Experimental Spectrum")

#select which reference spectra to use in the linear model
reference_spectra <- list(subset_reference_vectors[[6]],subset_reference_vectors[[7]], subset_reference_vectors[[9]], subset_reference_vectors[[12]])

#create data frame for references
ref <- as.data.frame(reference_spectra)
print(nrow(ref))
colnames(ref) <- c("CytbH", "CytbL", "Cytc1", "Cytc")
ref$Wavelength=wavelengths3
print(ncol(ref))

#creating data frame for spectra
subset_exp_df <- as.data.frame(subset_exp_vectors)

print(length(subset_exp_df))

formula <- response_variable ~ (I(ref$CytbH) * a) + 
  (I(ref$CytbL) * b) + 
  (I(ref$Cytc1) * c) + 
  (I(ref$Cytc) * d) + 
  (I(ref$Wavelength) * e + f)

initial_values <- list(a = 0.02, b = 0.03, c = 0.02, d=0.025, e = -0.02, f = 0)

#lower_bounds <- list(a = 0, b = 0, c = 0, d = 0, e = 0, f = 0, g = -Inf, h = -Inf)

model_results <- list()

#fitting over time
for (i in 1:(ncol(subset_exp_df))) {  
  #extract the response variable (experimental spectrum) for the current column
  response_variable <- subset_exp_df[[i]]
  
  #fit the nlsLM model for the current spectrum
  nls_model <- nlsLM(formula, 
                     data = cbind(response_variable, ref),
                     #lower=lower_bounds,
                     start = initial_values)
  
  #store the model results
  model_results[[colnames(subset_exp_df)[i]]] <- summary(nls_model)
}

print(model_results[1])

#the number of spectra after averaging should match the number of best fits
print(length(subset_exp_vectors))
print(length(model_results))

abs_bH <- list()  #initialize a list to store the fit results
abs_bL <- list()
abs_c1 <- list()
abs_c <- list()

#iterate over each model in model_results
for (i in seq_along(model_results)) {
  summary_obj <- model_results[[i]]
  coef_i <- coef(summary_obj)
  param_i <- coef_i[, 1]  
  
  #extract the peak index for each reference
  max_index_bH <- which.max(ref$CytbH)
  max_index_bL <- which.max(ref$CytbL)
  max_index_c1 <- which.max(ref$Cytc1)
  max_index_c <- which.max(ref$Cytc)
  peak_abs_bH <- ref$CytbH[max_index_bH]
  peak_abs_bL <- ref$CytbL[max_index_bL]
  peak_abs_c1 <- ref$Cytc1[max_index_c1]
  peak_abs_c <- ref$Cytc[max_index_c]
  
  #calculate absorbance by multiplying the reference spectra by all instances of 'a' coefficient
  abs_bH_i <- peak_abs_bH * param_i['a']  #assuming 'a' is the coefficient of interest
  abs_bL_i <- peak_abs_bL * param_i['b']
  abs_c1_i <- peak_abs_c1 * param_i['c']
  abs_c_i <- peak_abs_c * param_i['d']
  
  #store the result in the list
  abs_bH[[i]] <- abs_bH_i
  abs_bL[[i]] <- abs_bL_i
  abs_c1[[i]] <- abs_c1_i
  abs_c[[i]] <- abs_c_i
  
}

print(length(abs_bH))
num_vectors=(length(abs_bL))
time2 <- seq(0, num_vectors-1)
print(length(time2))
plot(time2, abs_bH, type = "l", xlab = "Time (sec)", ylab = "OD", main = paste("Cytochrome bH over time"))

abs_df <- data.frame(Time = time2, abs_bH = unlist(abs_bH),
                     abs_bL = unlist(abs_bL), abs_c1 = unlist(abs_c1),
                     abs_c = unlist(abs_c))

ggplot(abs_df, aes(x = Time)) +
  geom_line(aes(y = abs_bH, color = "Cytochrome bH")) +
  geom_line(aes(y = abs_bL, color = "Cytochrome bL")) +
  geom_line(aes(y = abs_c1, color = "Cytochrome c1")) +
  geom_line(aes(y = abs_c, color = "Cytochrome c")) +
  labs(x = "Time (seconds)", y = "OD", color = "Variable") +
  theme_minimal()

#calculate fbL for each time point using Map
fbL <- Map(function(bL, bH) bL / (bL + bH), abs_bL, abs_bH)

# Convert fbL to a data frame
fbL_df <- data.frame(Time = time2, fbL = unlist(fbL))

ggplot(fbL_df, aes(x = Time, y = fbL)) +
  geom_line(color = "blue") +
  labs(x = "Time (sec)", y = "fbL", color = "Spectrum") +
  theme_minimal()

save_data <- data.frame(Time = time2, abs_bH = unlist(abs_bH),
                              abs_bL = unlist(abs_bL), abs_c1 = unlist(abs_c1),
                              abs_c = unlist(abs_c), fbL=unlist(fbL))
print(file_path)

#create a descriptive file header (e.g., date_condition)
file_header <- readline("Enter the file header: ")

file_name <- paste(file_header, "save_data.csv", sep = "_")

#save the data frame to a CSV file with the constructed file name
write.csv(save_data, file = file_name, row.names = FALSE)

