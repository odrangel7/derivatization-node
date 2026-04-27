# EXPECTED COMPOUNDS 2.0 - SpectraMinds  
#===============================================================================
# Name   : SpectraMinds Expected Compounds (Registered Node Version)
# Author : Oscar Daniel Rangel Huerta (SpectraMinds), Optimized by AI
# Version: v4.0
# ===============================================================================

needed_pkgs <- c("data.table", "RJSONIO", "dplyr", "stringr", "tidyr")
for (p in needed_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(RJSONIO)

#------------------#
# Helper Functions #
#------------------#
safe_fread <- function(path) {
  if (!file.exists(path)) {
    log_msg(sprintf("Error: File not found: %s", path))
    stop(sprintf("File not found: %s", path))
  }
  tryCatch(
    {
      dt <- fread(path)
      log_msg(sprintf("Read %s: %d rows, %d cols", basename(path), nrow(dt), ncol(dt)))
      return(dt)
    },
    error = function(e) {
      log_msg(sprintf("Read error %s: %s", path, e$message))
      stop(sprintf("Read error %s: %s", path, e$message))
    }
  )
}

safe_write <- function(df, path) {
  tryCatch(
    {
      write.table(df, file = path, sep = "\t", row.names = FALSE, quote = TRUE)
      log_msg(sprintf("Wrote %s", basename(path)))
    },
    error = function(e) {
      log_msg(sprintf("Write error %s: %s", path, e$message))
      stop(sprintf("Write error %s: %s", path, e$message))
    }
  )
}

insert_col_desc <- function(cd, newcol) append(cd, list(newcol))

#---------------------------#
# 1. Argument & JSON Parsing #
#---------------------------#
# Read arguments and parse JSON input
args <- commandArgs()
inputFile <- args[6]
#inputFile <- "D:/Metabolomics/CompoundDiscoverer_Scripts/registered_node/node_args.json" # use instead of args[6] for testing with a local file

# Read JSON
json_content <- readLines(inputFile, warn = FALSE)
json_content <- paste(json_content, collapse = "\n")
json_content <- gsub("\u00A0", "", json_content) # Clean non-breaking spaces
CD_json_in <- fromJSON(json_content, simplify = FALSE)

# Setup Output Directory (Use input directory to avoid permission issues)
outdir <- file.path("c:/temp", "cd_out") ##I use this because sometimes CD has issues with writing on its own directories
if (!dir.exists(outdir)) dir.create(outdir, showWarnings = FALSE)

# Setup Logging
log_file <- file.path(outdir, "script_log.txt")
log_msg <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, msg), file = log_file, append = TRUE)
  message(msg)
}

# Clear previous log
if (file.exists(log_file)) file.remove(log_file)
log_msg("Script started.")
log_msg(sprintf("Input file: %s", inputFile))

result_out <- file.path(outdir, "result_out.csv")
cmpdsperfile_out <- file.path(outdir, "CmpdsPerFile_out.csv")
pairs_out <- file.path(outdir, "Pairs_table_out.csv")
pairs_long_out <- file.path(outdir, "Pairs_long_out.csv")


#---------------------------#
# 2. Build Reaction Config  #
#---------------------------#
params <- CD_json_in$NodeParameters

##read node parameters 
tol <- as.numeric(CD_json_in$NodeParameters[["Tolerance"]])

# Construct config table from UI parameters
# We check Thiol, DiThiol, Custom 1, Custom 2
potential_reactions <- list(
  list(Label = "Thiol", Mass = as.numeric(CD_json_in$NodeParameters[["Thiol Mass"]])),
  list(Label = "Di_Thiol", Mass = as.numeric(CD_json_in$NodeParameters[["DiThiol Mass"]])),
  list(Label = "Custom 1", Mass = as.numeric(CD_json_in$NodeParameters[["Custom 1 Mass"]]))
)

# Convert to data.table and filter out disabled (Mass = 0) reactions
reaction_config <- rbindlist(potential_reactions)
reaction_config <- reaction_config[Mass != 0 & !is.na(Mass)]

message(sprintf("Running with Tolerance=%.4f", tol))
message("Active Reactions:")
print(reaction_config)

#---------------------------#
# 3. Data Import            #
#---------------------------#
Compounds <- safe_fread(CD_json_in$Tables[[1]]$DataFile)
CmpdsPerFile <- safe_fread(CD_json_in$Tables[[2]]$DataFile)
InpFile <- safe_fread(CD_json_in$Tables[[3]]$DataFile)

setDT(Compounds)
setDT(CmpdsPerFile)
setDT(InpFile)

#---------------------------#
# 4. Processing Logic       #
#---------------------------#

# Calculate Expected Masses in CmpdsPerFile
for (i in seq_len(nrow(reaction_config))) {
  rxn <- reaction_config[i]
  col_name <- paste0("Expected_", rxn$Label)
  CmpdsPerFile[, (col_name) := `Calc MW` + rxn$Mass]
}

# Identify Controls
ControlSampleIDs <- InpFile[grepl("Control", `Sample Type`, ignore.case = TRUE), `Study File ID`]
ControlCmpds <- CmpdsPerFile[`Study File ID` %in% ControlSampleIDs]

# Loop Reactions for Matching
for (i in seq_len(nrow(reaction_config))) {
  rxn <- reaction_config[i]

  col_expected <- paste0("Expected_", rxn$Label)
  col_is_prod <- rxn$Label
  col_is_prec <- paste0("Suspected_precursor_", rxn$Label)

# Initialize the flag columns to FALSE if they don't exist yet to prevent crashes
if (!col_is_prod %in% names(Compounds)) Compounds[, (col_is_prod) := FALSE]
if (!col_is_prod %in% names(ControlCmpds)) ControlCmpds[, (col_is_prod) := FALSE]
  
  # A. Identify Products
  # Compounds matching Expected Mass from Controls
  ControlCmpds[, `:=`(Target_Min = get(col_expected) - tol, 
                      Target_Max = get(col_expected) + tol)]

  matched_idx <- Compounds[ControlCmpds,
    on = .(`Calc MW` >= Target_Min, `Calc MW` <= Target_Max),
    which = TRUE, nomatch = NULL
  ]

  Compounds[, (col_is_prod) := FALSE]
  Compounds[unique(matched_idx), (col_is_prod) := TRUE]

  # B. Identify Precursors
  # Get the products we just found
  prod_subset <- Compounds[get(col_is_prod) == TRUE]

  if (nrow(prod_subset) > 0) {
    # 1. Define bounds in the Product Subset
    prod_subset[, `:=`(MW_min = `Calc MW` - tol, MW_max = `Calc MW` + tol)]
    
    # We cannot use get(col_expected) inside 'on'. 
    # We create 'Temp_Join_Mass' in ControlCmpds to represent the dynamic column.
    ControlCmpds[, Temp_Join_Mass := get(col_expected)]
    
    # 3. Perform the Join
    matched_ctrl_idx <- ControlCmpds[prod_subset,
                                     on = .(Temp_Join_Mass >= MW_min, Temp_Join_Mass <= MW_max),
                                     which = TRUE, nomatch = NULL
    ]
    
    
    # 4. Clean up temp column immediately
    ControlCmpds[, Temp_Join_Mass := NULL]
    
    # 5. Update Control Flags
    ControlCmpds[, (col_is_prod) := FALSE]
    ControlCmpds[unique(matched_ctrl_idx), (col_is_prod) := TRUE]
    
  } else {
    ControlCmpds[, (col_is_prod) := FALSE]
  }
    
  # ---------------------------------------------------------
  # C. Map Valid Precursors back to Compounds
  # ---------------------------------------------------------
  
  # LOGIC FIX: Define valid_precursors AFTER the flags in ControlCmpds have been updated above
  valid_precursors <- ControlCmpds[get(col_is_prod) == TRUE]
  
  # Flag in CmpdsPerFile
  # Note: ensure 'Compounds per File ID' exists in both
  if(nrow(valid_precursors) > 0) {
    ids_to_flag <- valid_precursors$`Compounds per File ID`
    CmpdsPerFile[, (col_is_prec) := `Compounds per File ID` %in% ids_to_flag]
  } else {
    CmpdsPerFile[, (col_is_prec) := FALSE]
  }
  
  # Flag in Compounds (Precursor Flag)
  if (nrow(valid_precursors) > 0) {
    valid_precursors[, `:=`(MW_min = `Calc MW` - tol, MW_max = `Calc MW` + tol)]
    
    matched_prec_idx <- Compounds[valid_precursors,
                                  on = .(`Calc MW` >= MW_min, `Calc MW` <= MW_max),
                                  which = TRUE, nomatch = NULL
    ]
    Compounds[, (col_is_prec) := FALSE]
    Compounds[unique(matched_prec_idx), (col_is_prec) := TRUE]
  } else {
    Compounds[, (col_is_prec) := FALSE]
  }
  
  # Cleanup temp col
  rm_cols <- grep("Target_|MW_|Exp.*_min|Exp.*_max", names(ControlCmpds), value = TRUE)
  if (length(rm_cols) > 0) ControlCmpds[, (rm_cols) := NULL]
}
  
  
#---------------------------#
# 5. Pairs Logic            #
#---------------------------#
Compounds[, Pairs := ""]
pair_matches_list <- list()

for (i in seq_len(nrow(reaction_config))) {
  rxn <- reaction_config[i]
  col_prod <- rxn$Label
  col_prec <- paste0("Suspected_precursor_", rxn$Label)

  # Skip if cols missing (safety)
  if (!col_prod %in% names(Compounds) || !col_prec %in% names(Compounds)) next

  precursors <- Compounds[get(col_prec) == TRUE, .(ID = `Compounds ID`, MW = `Calc MW`)]
  products <- Compounds[get(col_prod) == TRUE, .(ID = `Compounds ID`, MW = `Calc MW`)]

  if (nrow(precursors) == 0 || nrow(products) == 0) next

  # Match: Product = Precursor + Mass
  precursors[, `:=`(Min = MW + rxn$Mass - tol, Max = MW + rxn$Mass + tol)]
  matches <- products[precursors, on = .(MW >= Min, MW <= Max), .(Prod_ID = x.ID, Prec_ID = i.ID), nomatch = NULL]

  matches <- matches[Prod_ID != Prec_ID]
  if (nrow(matches) > 0) {
    matches[, PairLabel := paste0(rxn$Label, "_", 1:.N)]
    pair_matches_list[[rxn$Label]] <- matches
  }
}

if (length(pair_matches_list) > 0) {
  all_pairs <- rbindlist(pair_matches_list)
  stack <- rbind(
    all_pairs[, .(ID = Prec_ID, Label = PairLabel)],
    all_pairs[, .(ID = Prod_ID, Label = PairLabel)]
  )
  labels_per_id <- stack[, .(NewPairs = paste(unique(Label), collapse = ", ")), by = ID]
  Compounds <- merge(Compounds, labels_per_id, by.x = "Compounds ID", by.y = "ID", all.x = TRUE)
  Compounds[!is.na(NewPairs), Pairs := NewPairs]
  Compounds[, NewPairs := NULL]
}

#---------------------------#
# 6. Save & Update JSON     #
#---------------------------#
safe_write(Compounds, result_out)
safe_write(CmpdsPerFile, cmpdsperfile_out)

# Construct Response JSON
node_response <- list(
  CurrentWorkflowID = CD_json_in$CurrentWorkflowID,
  ExpectedResponsePath = CD_json_in$ExpectedResponsePath,
  ResultFilePath = CD_json_in$ResultFilePath,
  NodeParameters = CD_json_in$NodeParameters,
  Version = CD_json_in$Version,
  Tables = CD_json_in$Tables
)

# -- Update Compounds (Table 1) --
cd1 <- node_response$Tables[[1]]$ColumnDescriptions
for (i in seq_len(nrow(reaction_config))) {
  rxn <- reaction_config[i]
  # Boolean Product
  cd1 <- insert_col_desc(cd1, list(
    ColumnName = rxn$Label, ID = "", DataType = "Boolean",
    Options = list(PositionAfter = "Calc MW", SpecialCellRenderer = "330D9522-50CC-41B7-9D45-5E5D8F708103")
  ))
  # Boolean Precursor
  cd1 <- insert_col_desc(cd1, list(
    ColumnName = paste0("Suspected_precursor_", rxn$Label), ID = "", DataType = "Boolean",
    Options = list(PositionAfter = "Calc MW", SpecialCellRenderer = "4B8956F1-3A54-423D-8654-60676825230F")
  ))
}
cd1 <- insert_col_desc(cd1, list(ColumnName = "Pairs", ID = "", DataType = "String", Options = list(PositionAfter = "Calc MW")))
node_response$Tables[[1]]$ColumnDescriptions <- cd1
node_response$Tables[[1]]$DataFile <- result_out

# -- Update CmpdsPerFile (Table 2) --
cd2 <- node_response$Tables[[2]]$ColumnDescriptions
for (i in seq_len(nrow(reaction_config))) {
  rxn <- reaction_config[i]
  cd2 <- insert_col_desc(cd2, list(
    ColumnName = paste0("Suspected_precursor_", rxn$Label), ID = "", DataType = "Boolean",
    Options = list(SpecialCellRenderer = "4B8956F1-3A54-423D-8654-60676825230F")
  ))
}
node_response$Tables[[2]]$ColumnDescriptions <- cd2
node_response$Tables[[2]]$DataFile <- cmpdsperfile_out


#---------------------------#
##We remove the tables 3, 4 and 5 to avoid problems when importing the data into the results table.
node_response$Tables <- node_response$Tables[1:2]

# Define 'ExpectedResponsePath' location.
json_out_file <- node_response$ExpectedResponsePath

# Write Response
response_json <- toJSON(node_response, indent = 1)
response_json <- gsub("\\[\\s*\\]", "{}", response_json)
jsonfileconn <- file(json_out_file) # jsonOutFile is the full path of the node_response.json file.  This needs to be in the same folder as node_args.json
writeLines(response_json, jsonfileconn)
close(jsonfileconn)


###save an image in the same directory to check on R
save.image(file = paste0(
  dirname(
    sub(".*=", "", args[4])
  ), "/CD node Rimage.dat"
))

message("Node execution complete.")
