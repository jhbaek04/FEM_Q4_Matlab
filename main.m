% Q4 FEM for plane-strain Linear Elasticity
% SE276C/MAE232C, Spring 2023, JS Chen
% Written by Jon Baek (job011@eng.ucsd.edu, jhbaek04@gmail.com)
clear; close all

%% set up the model
Model = model_setup;  % Set up your model in "model_setup.m" to be read here.

%% discretize
Mesh = sub_discretization ( Model );

%% detect boundary nodes
BC = sub_get_boundary ( Model, Mesh );

%% assembly
[ K , f ]  =  sub_assembly ( Model , Mesh , BC );

%% solve K * d = f
[ d ] = sub_solution ( K , f , BC );

%% postprocess
sub_postprocess ( Model , Mesh , d );
