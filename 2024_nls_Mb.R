install.packages("minpack.lm")
install.packages("ggplot2")
library(minpack.lm)
library(ggplot2)

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
  geom_line(aes(y = reference_vectors[[3]], color = "MbO")) +
  labs(x = "Wavelength", y = "OD", color = "References") +
  scale_color_manual(values = c("Cytochrome bH" = "blue", "MbO" = "red")) +
  theme_minimal()

#establish the bandwidth used for analyses
start_wavelength <- 535
end_wavelength <- 585

#plot absorbance over time
num_vectors <- length(abs_vectors)

#select the wavelength you want to evaluate over time
target_wavelength <- 580
target_index <- which.min(abs(wavelengths2 - target_wavelength))
print(target_index)

#extract intensities at the target wavelength from all spectra
intensities <- sapply(abs_vectors, function(vector) {
  #return the intensity at the target index
  vector[target_index]
})

#create a sequence of time points (assuming 1 spectrum per sec)
time <- seq(0, num_vectors-1)

#plot intensities as a function of time
plot(time, intensities, type = "l", xlab = "Time (seconds)", ylab = "Intensity", main = paste("Intensity at x over Time"))

intensity_time <- data.frame(time, intensities)
#write.csv(intensity_time, file = "time.csv", row.names = FALSE)

#find the indices corresponding to the start and end wavelengths
start_index <- which.min(abs(wavelengths2 - start_wavelength))
end_index <- which.min(abs(wavelengths2 - end_wavelength))
print(start_index)
print(end_index)

#manually enter the sample range you want to evaluate
vector_range = 3470:3520

#subset abs_vectors based on the user input
selected_vectors <- abs_vectors[vector_range]

#averages the spectra selected
exp_spectrum <- rowMeans(do.call(cbind, abs_vectors[vector_range]))
plot(wavelengths2, exp_spectrum, type = "l", xlab = "Wavelength", ylab = "OD", main = "Experimental Spectrum")

#crop the experimental spectrum and wavelength vector to the desired bandwidth
response_variable <- exp_spectrum[start_index:end_index]
wavelengths3=wavelengths2[start_index:end_index]

#extracts the peak wavelength of the experimental spectrum
max_index_spec <- which.max(response_variable)
peak_wav_spec=wavelengths3[max_index_spec]
print(peak_wav_spec)

#plot the cropped experimental spectrum against wavelength
plot(wavelengths3, response_variable, type = "l", xlab = "Wavelength (nm)", ylab = "OD", main = "Experimental Spectrum")

#crop the reference spectra to the desired bandwidth
subset_reference_vectors <- list()

for (i in seq_along(reference_vectors)) {
  subset_vector <- reference_vectors[[i]][start_index:end_index]
  subset_reference_vectors[[i]] <- subset_vector
}

plot(wavelengths3, subset_reference_vectors[[12]], type = "l", xlab = "Wavelength", ylab = "OD", main = "Cytochrome c")

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

max_index_MbD <- which.max(subset_reference_vectors[[2]])
peak_wav_MbD=wavelengths3[max_index_MbD]
print(peak_wav_MbD)

max_index_MbO <- which.max(subset_reference_vectors[[3]])
peak_wav_MbO=wavelengths3[max_index_MbO]
print(peak_wav_MbO)

#select which reference spectra to use in the linear model
#note, the original reference file has many unused columns
reference_spectra <- list(subset_reference_vectors[[2]], subset_reference_vectors[[3]], subset_reference_vectors[[6]],subset_reference_vectors[[7]], subset_reference_vectors[[9]], subset_reference_vectors[[12]])

#create new data frame
data <- as.data.frame(reference_spectra)
colnames(data) <- c("Column1", "Column2", "Column3", "Column4", "Column5", "Column6")
print(ncol(data))
data$ResponseVar=response_variable
data$Wavelength=wavelengths3
print(ncol(data))

formula <- response_variable ~ (Column1 * a) + (Column2 * b) + (Column3 * c) + (Column4 * d)+(Column5 * e)+(Column6 * f)+(Wavelength * g + h)
initial_values <- list(a=0.4, b=0.05, c = 0.02, d = 0.03, e = 0.02, f=0.025, g = -0.02, h = 0)

nls_model <- nlsLM(formula, data = data, start = initial_values)
summary(nls_model)
bestfit=predict(nls_model)

all_coef = coef(nls_model)
ref_coef=all_coef[1:(ncol(data))]
print(ref_coef)

fit_r1=data$Column1*ref_coef[1]
fit_r2=data$Column2*ref_coef[2]
plot(wavelengths3, fit_r2, type = "l", xlab = "Wavelength", ylab = "Intensity", main = "MbO")
fit_r3=data$Column3*ref_coef[3]
fit_r4=data$Column4*ref_coef[4]
fit_r5=data$Column5*ref_coef[5]
fit_r6=data$Column6*ref_coef[6]
line_fit=data$Wavelength*ref_coef[7]+ref_coef[8]

#calculate the absorbance of each cytochrome by multiplying the fit coefficient x the original reference at the peak wavelength
abs_bH<- fit_r3[max_index_bH]
abs_bL<- fit_r4[max_index_bL]
abs_c<- fit_r6[max_index_c]
abs_MbD<- fit_r1[max_index_MbD]
abs_MbO<- fit_r2[max_index_MbO]

#calculate fbL
fbL=abs_bL/(abs_bH+ abs_bL)

print(abs_bH)
print(abs_bL)
print(abs_c)
print(abs_MbD)
print(abs_MbO)
print(fbL)
#print(vector_range)

#calculate sum squared residuals
residuals <- residuals(nls_model)
squared_residuals <- residuals^2
sum_squared_residuals <- sum(squared_residuals)
print(sum_squared_residuals)

#correct spectrum and best fit
flatspec=response_variable-fit_r1-fit_r2-line_fit
flatfit=bestfit-line_fit
flatterfit=flatfit-fit_r1-fit_r2
plot(wavelengths3, flatspec, type = "l", xlab = "Wavelength", ylab = "OD", main = "Experimental Spectrum")


#plot residuals against wavelengths
plot(wavelengths3, residuals, type = "l", xlab = "Wavelength", ylab = "Residuals", main = "Residuals Plot")


#create data frame with experimental spectrum, best fit, and more
plot_data <- data.frame(Wavelength = wavelengths3, Experimental = response_variable, BestFit = bestfit, MbO=fit_r1, MbD=fit_r2, CytbH = fit_r3, CytbL=fit_r4, Cytc1=fit_r5, Cytc=fit_r6, Line=line_fit, Residual=residuals, FlatFit = flatfit, FlatterFit=flatterfit, FlatSpec = flatspec)

#plot the raw best fit
ggplot(plot_data, aes(x = Wavelength)) +
  geom_line(aes(y = Experimental), color = "blue") +
  geom_line(aes(y = BestFit), color = "red") +
  geom_line(aes(y = MbO), color = "black") +
  geom_line(aes(y = MbD), color = "grey") +
  geom_line(aes(y = CytbH), color = "pink") +
  geom_line(aes(y = CytbL), color = "purple") +
  geom_line(aes(y = Cytc1), color = "green") +
  geom_line(aes(y = Cytc), color = "orange") +
  geom_line(aes(y = Line), color = "yellow") +
  geom_line(aes(y = Residual), color = "magenta") +
  labs(x = "Wavelength", y = "OD", title = "Experimental Spectrum vs Best Fit") +
  theme_minimal()

#plot flat spec and flat fit
ggplot(plot_data, aes(x = Wavelength)) +
  geom_line(aes(y = FlatterFit), color = "blue") +
  geom_line(aes(y = flatspec), color = "red") +
  geom_line(aes(y = CytbH), color = "pink") +
  geom_line(aes(y = CytbL), color = "purple") +
  geom_line(aes(y = Cytc1), color = "green") +
  geom_line(aes(y = Cytc), color = "orange") +
  geom_line(aes(y = Residual), color = "magenta") +
  labs(x = "Wavelength", y = "OD", title = "Experimental Spectrum vs Best Fit") +
  theme_minimal()

print(file_path)

v#create a descriptive file header (e.g., date_condition)
file_header <- readline("Enter the file header: ")

file_name <- paste(file_header, "plot_data.csv", sep = "_")

#save the data frame to a CSV file with the constructed file name
write.csv(plot_data, file = file_name, row.names = FALSE)

