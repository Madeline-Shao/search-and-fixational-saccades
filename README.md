# search-and-fixational-saccades
Analyzes and plots eye tracker data, detects microsaccades


MicrosaccadeDetection
Filters data using filters function. Uses a velocity and acceleration threshold to mark microsaccades of each trial on a graph, then allows user to adjust it through adding and deleting microsaccades (through function edit_saccades). 

filters
Called by MicrosaccadeDetection, filters the data using either median, gaussian, or smoothing, then followed by sgolay.

edit_saccades
Called by MicrosaccadeDetection. User adds and deletes micro saccades.

FileProcessingCode
Processes a file into a structure containing fields for numerical data (image number, reaction time, condition, correctness, fixation duration, center coordinates, displacement, image/target/distractor coordinates, click coordinates, jitter) and alphabetical data (name, eye tested, scene, target)

TargetDetails
Finds height, width, and area in pixels and degrees of an image given coordinates of the upper left and lower right corners

find_stats
Finds the subject ID, percent correct, number incorrect, and average reaction time

data_stats
Uses find_stats to find the stats (subject ID, percent correct, number incorrect, average reaction time) of the data for each eye condition and saves it

eye_tracking_parse
Parses eye tracking data and removes blinks

training_data
Plots saccades in data using function find_saccades

find_saccades 
Called by training_data, plots data and allows user to add saccades

