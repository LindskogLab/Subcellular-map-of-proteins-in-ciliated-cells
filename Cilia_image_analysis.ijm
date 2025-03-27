// Title: 
// Analysis: Analyze overlap of candidate protein with panel markers in motile ciliated cells
// Author: Filippa Bertilsson
// Script written and run on: MacBook Pro M1 and M4
// Date: 2025-01-20

//mSet batch mode to let it run in the background
setBatchMode(true);

// Tidy up screen (if needed)
if (nImages>0) run("Close All"); // if there are 1 or more images - close all
print("\\Clear"); // empty log window
roiManager("reset"); // empty ROI manager
run("Options...", "iterations=1 count=1 black edm=32-bit"); // set Binary Options


// ANALYSIS
// First, set directory where you have your images
// Here there was one parent directory, with subfolers for each candidate protein
// Parent directory containing all protein folders
parentDir = "/Your/Directory/Here/";

// Get all subfolders (protein folders) in the parent directory
folders = getSubfolders(parentDir);


// Loop

// Main loop: This loop runs through each subfolder (protein folder), processes them, fins the .tif files, opens and analyzes them
// Loop through each subfolder
for (i = 0; i < folders.length; i++) {
    Folder = folders[i]; // Get path for a subfolder
    print("Processing folder: " + Folder);

    FileListArray = getFileList(Folder); // Get list of files (.tif) in folder
    Array.print(FileListArray);

    // Inner loop: This loop processes each .tif image
    for (w = 0; w < FileListArray.length; w++) {
        picname = FileListArray[w];
        print(picname);
        endsWithTif = endsWith(picname, ".tif"); // Checks that the file ends with ".tif", if so, it continues with processing 
        
        // If the file was a .tif, it opems it and extract info about it 
        if (endsWithTif == 1) {
            openPic = Folder + picname;
            options = "open=" + openPic + " autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_1";
            run("Bio-Formats Importer", options); // Opens the file with Bio-Formats

            // Get image filename without extension (.tif)
            image = getInfo("image.filename");
            dotIndex = indexOf(image, ".tif" );
            imageWithoutExtension = substring(image, 0, dotIndex); // This line removes ".tif"

            // Extract prefix 
            // This part extracts the Slide ID (e.g., "20465761") from the file name, this is how we kept trach of which protein was stained with which slide
            underscoreIndex = indexOf(imageWithoutExtension, "_");
            if (underscoreIndex > -1) {
                imagePrefix = substring(imageWithoutExtension, 0, underscoreIndex);
            } else {
                imagePrefix = "UnknownPrefix";
            }


            // Rename files based on the specified patterns
            // Each core on the TMA had a coordinate (e.g., [1,1,A], telling us which tissue it was
            // This part of the loop identifies which tissue it is for each image being analyzed, based on our coordinate system
            // Once identified, it renames the file using the tissue name, and the core ID (e.g., 1A)
            if (indexOf(imageWithoutExtension, "[1,1,A]") > -1) {
                newFileName = "Fallopian_1A";
            } else if (indexOf(imageWithoutExtension, "[1,1,B]") > -1) {
                newFileName = "Fallopian_1B";
            } else if (indexOf(imageWithoutExtension, "[1,2,A]") > -1) {
                newFileName = "Fallopian_2A";
            } else if (indexOf(imageWithoutExtension, "[1,2,B]") > -1) {
                newFileName = "Fallopian_2B";
            } else if (indexOf(imageWithoutExtension, "[1,3,A]") > -1) {
                newFileName = "Fallopian_3A";
            } else if (indexOf(imageWithoutExtension, "[1,3,B]") > -1) {
                newFileName = "Fallopian_3B";
            } else if (indexOf(imageWithoutExtension, "[1,1,C]") > -1) {
                newFileName = "Cervix_1C";
            } else if (indexOf(imageWithoutExtension, "[1,1,D]") > -1) {
                newFileName = "Cervix_1D";
            } else if (indexOf(imageWithoutExtension, "[1,2,C]") > -1) {
                newFileName = "Cervix_2C";
            } else if (indexOf(imageWithoutExtension, "[1,2,D]") > -1) {
                newFileName = "Cervix_2D";
            } else if (indexOf(imageWithoutExtension, "[1,3,C]") > -1) {
                newFileName = "Cervix_3C";
            } else if (indexOf(imageWithoutExtension, "[1,3,D]") > -1) {
                newFileName = "Cervix_3D";
            } else if (indexOf(imageWithoutExtension, "[1,4,A]") > -1) {
                newFileName = "Endometrium_4A";
            } else if (indexOf(imageWithoutExtension, "[1,4,B]") > -1) {
                newFileName = "Endometrium_4B";
            } else if (indexOf(imageWithoutExtension, "[1,5,A]") > -1) {
                newFileName = "Endometrium_5A";
            } else if (indexOf(imageWithoutExtension, "[1,5,B]") > -1) {
                newFileName = "Endometrium_5B";
            } else if (indexOf(imageWithoutExtension, "[1,6,A]") > -1) {
                newFileName = "Endometrium_6A";
            } else if (indexOf(imageWithoutExtension, "[1,6,B]") > -1) {
                newFileName = "Endometrium_6B";
            } else if (indexOf(imageWithoutExtension, "[1,4,C]") > -1) {
                newFileName = "Nasopharynx_4C";
            } else if (indexOf(imageWithoutExtension, "[1,4,D]") > -1) {
                newFileName = "Nasopharynx_4D";
            } else if (indexOf(imageWithoutExtension, "[1,5,C]") > -1) {
                newFileName = "Nasopharynx_5C";
            } else if (indexOf(imageWithoutExtension, "[1,5,D]") > -1) {
                newFileName = "Nasopharynx_5D";
            } else if (indexOf(imageWithoutExtension, "[1,6,C]") > -1) {
                newFileName = "Nasopharynx_6C";
            } else if (indexOf(imageWithoutExtension, "[1,6,D]") > -1) {
                newFileName = "Nasopharynx_6D";
            } else if (indexOf(imageWithoutExtension, "[1,7,C]") > -1) {
                newFileName = "Bronchus_7C";
            } else if (indexOf(imageWithoutExtension, "[1,7,D]") > -1) {
                newFileName = "Bronchus_7D";
            } else if (indexOf(imageWithoutExtension, "[1,8,C]") > -1) {
                newFileName = "Bronchus_8C";
            } else if (indexOf(imageWithoutExtension, "[1,8,D]") > -1) {
                newFileName = "Bronchus_8D";
            } else if (indexOf(imageWithoutExtension, "[1,9,C]") > -1) {
                newFileName = "Bronchus_9C";
            } else if (indexOf(imageWithoutExtension, "[1,9,D]") > -1) {
                newFileName = "Bronchus_9D";
            }

            print("New file name: " + newFileName);

            run("8-bit"); // Converts image to 8-bit 
            run("Stack to Images"); // Splits the stack into separate images

            // Reset contrast and brightness using function below
            processChannels(imageWithoutExtension);

            // Process ROI and measure using function below
            processROIs(imageWithoutExtension);
			
			// Save results, here the file names were named accordingly: "Tissue_CoreID_SlideID.csv"
            saveAs("Results", Folder + newFileName + "_" + imagePrefix + ".csv");
            
			// Clean up, prep for next folder to analyze
            roiManager("reset"); // empty ROI manager
            run("Options...", "iterations=1 count=1 black edm=32-bit"); // set Binary Options
            if (nImages > 0) run("Close All"); // if there are 1 or more images - close all
            if (isOpen("Results")) { // close results window, if open
                selectWindow("Results");
                close("Results");
            }
        }
    }
}

// Function to process image channels
// This function resets brightness and contrast and closes the images we do not need
function processChannels(imageName) {
    channels = newArray("-0007", "-0006", "-0005", "-0004", "-0003", "-0002");
    for (c = 0; c < channels.length; c++) {
        selectWindow(imageName + channels[c]);
        resetMinAndMax();
    }
    close(imageName + "-0001"); // DAPI
    close(imageName + "-0008"); // Autofluoresence 
    run("Tile");
}

// Function to process ROIs
// This function uses Multi OtsuThreshold to threshold the images of the panel markers ("Cytoplasm", "BB", "TZ", "AX", "Nucleus")
// If it is uncable to threshold an image, it skips to the next stack (tissue core)
// Using the thresholding, it generate ROIs which are then used to measure the overlap between each panel marker and the candidate protein
function processROIs(imageName) {
    channelNames = newArray("Cytoplasm", "RL", "TZ", "CL", "Nucleus"); // Labels for panel markers
    channelSuffixes = newArray("-0004", "-0006", "-0002", "-0005", "-0007"); // Panel marker identifiers 

    for (c = 0; c < channelNames.length; c++) {
        selectImage(imageName + channelSuffixes[c]);
        run("Multi OtsuThreshold", "numlevels=3"); // Thresholding / Segmentation

        // Close unused regions and keep only Region 2, the level with the highest pixel values, thus presumably most precise target signal
        selectImage("Region 0"); close();
        selectImage("Region 1"); close();
        selectImage("Region 2"); rename(channelNames[c]);
        selectImage("Histogram"); close();

        // Generate ROI based on segmentation (skip if it fails)
        setThreshold(1, 255);
        run("Create Selection");

        // Check if a selection was successfully created
        if (selectionType() == -1) { // No selection
            print("No ROI created for " + channelNames[c] + ". Skipping this file.");
            roiManager("reset"); // Clear any existing ROIs for the file
            if (nImages > 0) run("Close All"); // Close all images for the file
            return; // Skip the current file (stack) and move on to the next
        }
		
		// Add ROI to ROI manager
        roiManager("Add");
        roiManager("Select", roiManager("count") - 1); // Select the last-added ROI
        roiManager("Rename", channelNames[c]); // Renames each ROI to its panel marker label
    }

    // Use all ROIs to measure overlap with candidate protein image (only if ROI creation succeeded)
    selectImage(imageName + "-0003"); // Image with the candidate protein image
    for (r = 0; r < channelNames.length; r++) {
        roiManager("Select", r);
        roiManager("Measure"); // Measures: Area of ROI, Mean pixel value, Min pixel value, and Max pixel value
    }
}


// Helper functions 
// Function to get the name of the folder
function getFolderName(folderPath) {
    pathParts = split(folderPath, "/");
    return pathParts[pathParts.length - 2]; // Get folder name from path.
}

// Function to get subfolders of a directory
function getSubfolders(directory) {
    fileList = getFileList(directory);
    subfolders = newArray();
    for (i = 0; i < fileList.length; i++) {
        if (File.isDirectory(directory + fileList[i])) {
            subfolders = Array.concat(subfolders, newArray(directory + fileList[i] + "/"));
        }
    }
    return subfolders;
}


// Save the log manually if you want/need it:) 

