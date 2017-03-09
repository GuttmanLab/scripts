// Ask for the directory
dir = getDirectory("Choose a directory");
files = getFileList(dir);

// Dialog box options
projections = newArray("Max Intensity", "Average Intensity", "Min Intensity", "Sum Slices", "Standard Deviation", "Median");
colors = newArray("Blue", "Green", "Red", "Cyan", "Yellow", "Magenta");

// Create the dialog box
Dialog.create("User-defined variables");
Dialog.addNumber("Number of channels:", 3);
Dialog.addChoice("Z-projection: ", projections, "Max Intensity");
Dialog.addChoice("Channel 1 color:", colors, "Blue");
Dialog.addChoice("Channel 2 color:", colors, "Green");
Dialog.addChoice("Channel 3 color:", colors, "Red");
Dialog.addChoice("Channel 4 color:", colors, "Cyan");
Dialog.addChoice("Channel 5 color:", colors, "Yellow");
Dialog.addChoice("Channel 6 color:", colors, "Magenta");
Dialog.show();

// Read the dialog box
numChannels = Dialog.getNumber();
projection = Dialog.getChoice();
color1 = Dialog.getChoice();
color2 = Dialog.getChoice();
color3 = Dialog.getChoice();
color4 = Dialog.getChoice();
color5 = Dialog.getChoice();
color6 = Dialog.getChoice();

// ImageJ uses these three-letter codes to refer to different project schemes
// They appear in window names, and so are necessary to select windows.
if (projection == "Max Intensity") {
	proj = "MAX";
} else if (projection == "Min Intensity") {
	proj = "MIN";
} else if (projection == "Average Intensity") {
	proj = "AVG";
} else if (projection == "Median") {
	proj = "MED";
} else if (projection == "Standard Deviation") {
	proj = "STD";
} else {
	proj = "SUM";
}

setBatchMode(true);

// For each file...
for (f = 0; f < files.length; f++) {
    file = files[f];
    print("Processing " + file);

    // If the file ends in ".tif"
    if (endsWith(file, ".tif")) {
        baseNameEnd = indexOf(file, ".tif");
        baseName = substring(file, 0, baseNameEnd);

        // Open it
        run("Bio-Formats Importer", "open=" + dir + file + " autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");

        // Convert from a single stack (single channel) to a hyperstack (multiple channel)
        numSlicesPerChannel = nSlices / numChannels;       
        run("Stack to Hyperstack...", "order=xyzct channels=" + numChannels + " slices=" + numSlicesPerChannel + " frames=1 display=Color");

        // Set the channels of each channel
        colorScheme = newArray(color1, color2, color3, color4, color5, color6);
        for (i = 0; i < numChannels; i++) {
            Stack.setChannel(i + 1);
            run(colorScheme[i]);
        }

        // Z-project the image
        run("Z Project...", "projection=[" + projection + "]");
        selectWindow(proj + "_" + baseName + ".tif");
        saveAs("Tiff", dir + baseName + "_" + proj + ".tif");

        // Split the image into its individual channels and save each separately
        run("Split Channels");
        for (i = numChannels; i > 0; i--) {
            selectWindow("C" + i + "-" + baseName + "_" + proj + ".tif");
            run("8-bit");
            saveAs("Tiff", dir + baseName + "_" + proj + "_C" + i + ".tif");
            close();
        }

    close();
    }
}

setBatchMode(false);
