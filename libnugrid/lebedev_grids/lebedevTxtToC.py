import os
import glob
from math import pi
from math import sin
from math import cos

output_c_file = open("lebedevData.cpp", 'w')
output_c_file.write('#include "lebedevData.h"\n')
output_c_file.write("std::vector<std::vector<double> > getLebedevGridByNumPoints(int numPoints) {\n")
output_c_file.write("std::vector<double> xArray;\n")
output_c_file.write("std::vector<double> yArray;\n")
output_c_file.write("std::vector<double> zArray;\n")
output_c_file.write("std::vector<double> wArray;\n")
output_c_file.write("std::vector<std::vector<double> > returnArray;\n")
output_c_file.write("switch (numPoints) {\n")
f_list = glob.glob(os.getcwd() + "/*.txt")
f_list.sort()
for f in f_list:
	print f
	num_lines = sum(1 for line in open(f))
	output_c_file.write("\tcase {}: \n".format(num_lines))
	x_vector_string = ""
	y_vector_string = ""
	z_vector_string = ""
	w_vector_string = ""
	with open(f) as grid_file:
		for line in grid_file:
			point = line.split()
			theta = float(point[0]) / 180.0 * pi
			phi = float(point[1]) / 180.0 * pi
			x_vector_string += "{},".format(cos(theta) * sin(phi))
			y_vector_string += "{},".format(sin(theta) * sin(phi))
			z_vector_string += "{},".format(cos(phi))
			w_vector_string += "{},".format(float(point[2]) * 4 * pi) 
			
	output_c_file.write("xArray = {" + x_vector_string[:-1] + "};\n")
	output_c_file.write("yArray = {" + y_vector_string[:-1] + "};\n")
	output_c_file.write("zArray = {" + z_vector_string[:-1] + "};\n")
	output_c_file.write("wArray = {" + w_vector_string[:-1] + "};\n")
	output_c_file.write("returnArray = {xArray, yArray, zArray, wArray};\n")
	output_c_file.write("return returnArray;\n")
	
output_c_file.write("default:\n")
output_c_file.write("throw 'numer of points not in file';\n")
output_c_file.write("}\n")
output_c_file.write("}\n")
output_c_file.write("std::vector<std::vector<double> > getLebedevGridByNumPrecision(int precision) {\n")
output_c_file.write("std::vector<double> xArray;\n")
output_c_file.write("std::vector<double> yArray;\n")
output_c_file.write("std::vector<double> zArray;\n")
output_c_file.write("std::vector<double> wArray;\n")
output_c_file.write("std::vector<std::vector<double> > returnArray;\n")
output_c_file.write("switch (precision) {\n")
f_list = glob.glob(os.getcwd() + "/*.txt")
f_list.sort()
for f in f_list:
	print f
	num_lines = sum(1 for line in open(f))
	output_c_file.write("\tcase {}: \n".format(int(f[-7:-4])))
	x_vector_string = ""
	y_vector_string = ""
	z_vector_string = ""
	w_vector_string = ""
	with open(f) as grid_file:
		for line in grid_file:
			point = line.split()
			theta = float(point[0]) / 180.0 * pi
			phi = float(point[1]) / 180.0 * pi
			x_vector_string += "{},".format(cos(theta) * sin(phi))
			y_vector_string += "{},".format(sin(theta) * sin(phi))
			z_vector_string += "{},".format(cos(phi))
			w_vector_string += "{},".format(float(point[2]) * 4 * pi) 
			
	output_c_file.write("xArray = {" + x_vector_string[:-1] + "};\n")
	output_c_file.write("yArray = {" + y_vector_string[:-1] + "};\n")
	output_c_file.write("zArray = {" + z_vector_string[:-1] + "};\n")
	output_c_file.write("wArray = {" + w_vector_string[:-1] + "};\n")
	output_c_file.write("returnArray = {xArray, yArray, zArray, wArray};\n")
	output_c_file.write("return returnArray;\n")
output_c_file.write("default:\n")
output_c_file.write("throw 'numer of points not in file';\n")
output_c_file.write("}\n")
output_c_file.write("}\n")
