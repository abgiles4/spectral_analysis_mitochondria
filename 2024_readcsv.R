print("Choose the data you'd like to analyze:")
#the columns should be arranged as follows: wavelength (1), dark (2), IO (3), data (4...)
file_path <- file.choose()
print(file_path)
data <- read.csv(file_path)

#extract wavelengths and spectra vectors
wavelengths <- data[, 1]
dark=data[,2]
io=data[,3]
spectra_vectors <- lapply(data[, -(1:3)], as.vector)
num_spectra <- length(spectra_vectors)
plot(wavelengths, spectra_vectors[[9]], type = "l", xlab = "Wavelength (nm)", ylab = "Intensity", main = "Transmission")

#subtract dark from IO
io_corrected <- io - dark
plot(wavelengths,io_corrected,type = "l", xlab = "Wavelength (nm)", ylab = "Intensity", main = "Incident Light")

#subtract dark from each vector in spectra_vectors
spectra_corrected <- lapply(spectra_vectors, function(spectrum) spectrum - dark)
plot(wavelengths, spectra_corrected[[9]], type = "l", xlab = "Wavelength", ylab = "Intensity", main = "Transmission")

#calculating absorbance
log_ratio <- lapply(spectra_corrected, function(spectrum) log10(io_corrected / spectrum))

#plot the 9th absorbance spectrum against wavelength
plot(wavelengths, log_ratio[[9]], type = "l", xlab = "Wavelength (nm)", ylab = "OD", main = "Absorbance Spectrum")

#creating a new data.frame to store the data
log_ratio_df <- data.frame(wavelengths, log_ratio)

print(file_path)

#create a descriptive file header (e.g., date_condition)
file_header <- readline("Enter the file header: ")
file_name <- paste(file_header, "log_ratios.csv", sep = "_")

#save the data frame to a CSV file with the constructed file name
write.csv(log_ratio_df, file = paste0(file_name), row.names = FALSE)

