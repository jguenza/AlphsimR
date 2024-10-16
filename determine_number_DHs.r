



# Define the number of chromosomes and markers
num_chromosomes <- 7
markers_per_chromosome <- 20  # Adjusted for finer resolution

# Create marker names
markerNames <- paste0("M", 1:(num_chromosomes * markers_per_chromosome))

# Create chromosome assignments
chromosomes <- rep(1:num_chromosomes, each = markers_per_chromosome)

# Calculate positions based on recombination rate (3.2 cM per Mb)
# Assuming a physical distance of 1 Mb = 3.2 cM, markers are spaced ~0.3125 Mb apart
positions <- rep(seq(from = 0, to = 19 * 0.3125, by = 0.3125), times = num_chromosomes)

# Combine into a data frame
genMap <- data.frame(
  markerName = markerNames,
  chromosome = chromosomes,
  position = positions
)

# Display the genetic map
print(genMap)

# Adjust genotypes for 140 loci
geno = rbind(rep(2, 140), # Parent 1: All loci have genotype 2
             rep(0, 140)) # Parent 2: All loci have genotype 0

# Verify dimensions
print(dim(geno)) # Should print 2 rows by 140 columns

# Create pedigree with just IDs
# RP stands for recurrent parent
ped = data.frame(
  id = c("RP", "Donor"),
  mother = c(0, 0),
  father = c(0, 0)
)


# Create initial founder population
founderPop = importInbredGeno(geno = geno,
                              genMap = genMap,
                              ped = ped)
