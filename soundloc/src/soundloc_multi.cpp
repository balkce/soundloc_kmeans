/*
	Source coded by Caleb Rascon, 2010
	IIMAS, UNAM
	México
	
	Example that uses libmultisoundloc.a
*/

#include "multisoundloc.h"

#include <iostream>
#include <signal.h>
#include <string.h>
#include <stdio.h>

//ROS libs
#include "ros/ros.h"
#include "std_srvs/Empty.h"
#include <jsoncpp/json/json.h>
#include <json_msgs/JsonSrv.h>

//Glib, for threading (GThread)
#include <glib.h>

void finalize(int signum){
	std::cout << "\nSoundLoc: Caught signal " << signum << std::endl;
	std::cout << "SoundLoc: Exiting gracefully.\n";
	
	//Telling all threads that we are shutting down.
	
	soundloc_finalize();
	
	//gplot_stream.cmd("set terminal postscript");
	//gplot_stream.cmd("set output 'streams.eps'");
	//gplot_stream.cmd("replot");
	
	std::cout << "SoundLoc: Bye.\n";
	
	exit(EXIT_SUCCESS);
}

bool ros_reset_soundloc(json_msgs::JsonSrv::Request  &req, json_msgs::JsonSrv::Response &res){
	ROS_INFO("resetting soundloc");
	soundloc_clear();
	return true;
}

bool ros_get_sound_directions(json_msgs::JsonSrv::Request  &req, json_msgs::JsonSrv::Response &res){
	ROS_INFO("providing current sound directions");
	Json::FastWriter writer;

	Json::Value response;
	Json::Value json_sound_direction(Json::arrayValue);	
	if(sources.size() > 0){
		bool first_angle = TRUE;
		for(int i = 0; i <sources.size(); i++){
			Json::Value object_json;
			object_json["doa"] = sources[i].doa;
			json_sound_direction.append(object_json);
		}
    	response["doas"] = json_sound_direction;
	}else{
    	response = json_sound_direction;
	}

	res.json = writer.write(response);
	
	return true;
}

int main( int argc, char *argv[] )
{
	//if ( argc != 8 && argc != 9 && argc != 10) usage();
	
	//Assigning a CTRL+C handler
	signal(SIGINT, finalize);
	
	//ROS initialization
	std::cout << "SoundLoc: Starting up ROS connection...\n";
	
	ros::init(argc, argv, "soundloc");
	ros::NodeHandle ros_nh("~");
	ros::ServiceServer reset_soundloc = ros_nh.advertiseService("reset_soundloc", ros_reset_soundloc);
	ros::ServiceServer get_sound_directions = ros_nh.advertiseService("get_sound_directions", ros_get_sound_directions);
	
	// Obtaining parameters from ROS parameter server
	float distance_between_mics;
	if (ros_nh.getParam("distance_between_mics",distance_between_mics)){
		ROS_INFO("Distance between mics: %f",distance_between_mics);
	}else{
		distance_between_mics = 0.21;
		ROS_WARN("Distance between mics argument not found in ROS param server, using default value(%f).",distance_between_mics);
	}
	
	bool graph_out;
	if (ros_nh.getParam("graph_out",graph_out)){
		ROS_INFO("Graph out: %d",graph_out);
	}else{
		graph_out = true;
		ROS_WARN("Graph out argument not found in ROS param server, using default value(%d).",graph_out);
	}
	
	int max_number_sources;
	if (ros_nh.getParam("max_number_sources",max_number_sources)){
		ROS_INFO("Maximum number of sources: %d",max_number_sources);
	}else{
		max_number_sources = 3;
		ROS_WARN("Maximum number of sources argument not found in ROS param server, using default value(%d).",max_number_sources);
	}
	
	bool connect_ports;
	if (ros_nh.getParam("connect_ports",connect_ports)){
		ROS_INFO("Automatic port connection: %d",connect_ports);
	}else{
		connect_ports = false;
		ROS_WARN("Automatic port connection argument not found in ROS param server, using default value(%d).",connect_ports);
	}
	
	int gcc_style;
	if (ros_nh.getParam("gcc_style",gcc_style)){
		ROS_INFO("GCC style: %d",gcc_style);
	}else{
		gcc_style = 4;
		ROS_WARN("GCC style argument not found in ROS param server, using default value(%d).",gcc_style);
	}
	
	double gcc_th;
	if (ros_nh.getParam("gcc_th",gcc_th)){
		ROS_INFO("Correlation treshold: %f",gcc_th);
	}else{
		gcc_th = 100.0;
		ROS_WARN("Correlation treshold argument not found in ROS param server, using default value(%f).",gcc_th);
	}
	
	double redundancy_th;
	if (ros_nh.getParam("redundancy_th",redundancy_th)){
		ROS_INFO("DOA redundancy treshold: %f",redundancy_th);
	}else{
		redundancy_th = 20.0;
		ROS_WARN("DOA redundancy treshold argument not found in ROS param server, using default value(%f).",redundancy_th);
	}
	
	int dynamic_gcc_th;
	if (ros_nh.getParam("dynamic_gcc_th",dynamic_gcc_th)){
		ROS_INFO("Dynamic GCC: %d",dynamic_gcc_th);
	}else{
		dynamic_gcc_th = 1;
		ROS_WARN("Dynamic GCC argument not found in ROS param server, using default value(%d).",dynamic_gcc_th);
	}
	
	int moving_average;
	if (ros_nh.getParam("moving_average",moving_average)){
		ROS_INFO("Moving average: %d",moving_average);
	}else{
		moving_average = 1;
		ROS_WARN("Moving average argument not found in ROS param server, using default value(%d).",moving_average);
	}
	
	int moving_factor;
	if (ros_nh.getParam("moving_factor",moving_factor)){
		ROS_INFO("Moving factor: %d",moving_factor);
	}else{
		moving_factor = 1;
		ROS_WARN("Moving factor argument not found in ROS param server, using default value(%d).",moving_factor);
	}
	
	int memory_factor;
	if (ros_nh.getParam("memory_factor",memory_factor)){
		ROS_INFO("Memory factor: %d",memory_factor);
	}else{
		memory_factor = 5;
		ROS_WARN("Memory factor argument not found in ROS param server, using default value(%d).",memory_factor);
	}
	
	int kmeans_min_dist;
	if (ros_nh.getParam("kmeans_min_dist",kmeans_min_dist)){
		ROS_INFO("kmean minimum distance: %d",kmeans_min_dist);
	}else{
		kmeans_min_dist = 10;
		ROS_WARN("kmean minimum distance argument not found in ROS param server, using default value(%d).",kmeans_min_dist);
	}
	
	int max_plot_confidence;
	if (ros_nh.getParam("max_plot_confidence",max_plot_confidence)){
		ROS_INFO("Cluster size for maximum confidence: %d",max_plot_confidence);
	}else{
		max_plot_confidence = 4;
		ROS_WARN("Cluster size for maximum confidence argument not found in ROS param server, using default value(%d).",max_plot_confidence);
	}
	
	double noise_threshold;
	if (ros_nh.getParam("noise_threshold",noise_threshold)){
		ROS_INFO("Noise treshold: %f",noise_threshold);
	}else{
		noise_threshold = 0.001;
		ROS_WARN("Noise treshold argument not found in ROS param server, using default value(%f).",noise_threshold);
	}
	
	double noise_peak_change;
	if (ros_nh.getParam("noise_peak_change",noise_peak_change)){
		ROS_INFO("Noise peak change: %f",noise_peak_change);
	}else{
		noise_peak_change = 0.0015;
		ROS_WARN("Noise peak change argument not found in ROS param server, using default value(%f).",noise_peak_change);
	}
	
	bool verbose;
	if (ros_nh.getParam("verbose",verbose)){
		ROS_INFO("Verbose: %d",verbose);
	}else{
		verbose = false;
		ROS_WARN("Verbose argument not found in ROS param server, using default value(%d).",verbose);
	}
	
	std::cout << "\nSoundLoc: Stand by while everything is initialized.\n\n";
	
	soundloc_init(distance_between_mics, max_number_sources, graph_out, connect_ports, gcc_style, gcc_th, redundancy_th, dynamic_gcc_th, moving_average, moving_factor, memory_factor, kmeans_min_dist, max_plot_confidence, noise_threshold, noise_peak_change, verbose);
	
	ros::spin();
	return 0;
}
