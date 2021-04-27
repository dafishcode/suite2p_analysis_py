import numpy as np
import random
import matplotlib.pyplot as plt
import itertools

def create_empty_array(size):
    new_array = create_empty_list(size)

    for row in range(size):
        columns = create_empty_list(size)
        new_array[row] = columns

    return new_array

def create_zero_array(size):
    new_array = create_empty_list(size)

    for row in range(size):
        columns = create_zero_list(size)
        new_array[row] = columns

    new_array = np.array(new_array)
    return new_array


def create_irregular_zero_array(height,width):
    new_array = create_empty_list(height)

    for row in range(height):
        columns = create_zero_list(width)
        new_array[row] = columns

    new_array = np.array(new_array)
    return new_array


def create_empty_list(size):
    new_list = []
    for number in range(size):
        new_list.append([])

    return new_list

def create_zero_list(size):
    new_list = []
    for number in range(size):
        new_list.append(float(0))

    return new_list


def upper_triangle_to_symmetrical_matrix(original_matrix):

    matrix_size = len(original_matrix[0])
    new_matrix = create_empty_array(matrix_size)

    for osscilator_from in range(matrix_size):
        for osscilator_to in range(matrix_size):
            weight_1 = original_matrix[osscilator_from][osscilator_to]
            weight_2 = original_matrix[osscilator_to][osscilator_from]
            new_matrix[osscilator_from][osscilator_to] = max([weight_1,weight_2])

    new_matrix = np.array(new_matrix)
    return new_matrix


def symmetrical_to_upper_triangle(original_matrix):

    number_of_rows,number_of_columns = np.shape(original_matrix)
    new_matrix = create_zero_array(number_of_columns)

    for row in range(number_of_rows):
        for column in range(row, number_of_columns):
            value = original_matrix[row][column]
            new_matrix[row][column] = value

    return new_matrix


def normalise_connectivity_matrix(original_matrix):
    matrix_size, matrix_size2 = np.shape(original_matrix)
    normalised_connectivity_matrix = create_empty_array(matrix_size)
    maximum_weight = np.max(original_matrix)

    for oscillator_from in range(matrix_size):
        for oscillator_to in range(matrix_size):
            original_value = original_matrix[oscillator_from][oscillator_to]
            scaled_value = original_value / maximum_weight
            normalised_connectivity_matrix[oscillator_from][oscillator_to] = scaled_value

    normalised_connectivity_matrix = np.array(normalised_connectivity_matrix)

    return normalised_connectivity_matrix


def p_value_threshold_correlation_matrix(p_matrix,correlation_matrix,threshold):
    matrix_size,matrix_size_2 = np.shape(p_matrix)
    new_matrix = create_zero_array(matrix_size)

    for row in range(matrix_size):
        for column in range(matrix_size):

            if p_matrix[row][column] < threshold:
                new_matrix[row][column] = correlation_matrix[row][column]

    return new_matrix



def get_subgraph_list_connection_matrix(matrix,subnode_list):
    matrix_size = len(matrix[0])

    for row in range(matrix_size):
        for column in range(matrix_size):

            if row not in subnode_list or column not in subnode_list:
                list_size = len(matrix[row][column])
                matrix[row][column] = create_zero_list(list_size)
                matrix[column][row] = create_zero_list(list_size)

    return matrix



def get_subgraph_connection_matrix(matrix,subnode_list):
    matrix_size = len(matrix[0])

    for row in range(matrix_size):
        for column in range(matrix_size):

            if row not in subnode_list or column not in subnode_list:
                matrix[row][column] = 0
                matrix[column][row] = 0

    return matrix

def absolute_value_matrix(matrix):
    rows,columns = np.shape(matrix)

    for row in range(rows):
        for column in range(columns):
            original_value = matrix[row][column]
            new_value = np.abs(original_value)
            matrix[row][column] = new_value

    return matrix


def extract_submatrix(full_matrix,list_of_subnodes):

    number_of_subnodes = len(list_of_subnodes)
    submatrix = create_empty_array(number_of_subnodes)

    for subnode_index in range(number_of_subnodes):
        subnode = list_of_subnodes[subnode_index]

        for other_subnode_index in range(number_of_subnodes):
            other_subnode = list_of_subnodes[other_subnode_index]

            weight = full_matrix[subnode][other_subnode]
            submatrix[subnode_index][other_subnode_index] = weight

    submatrix = np.array(submatrix)
    return submatrix


def embed_submatrix(submatrix,full_matrix):

    submatrix_size = len(submatrix[0])

    for subnode1 in range(submatrix_size):
        for subnode2 in range(submatrix_size):
            weight = submatrix[subnode1][subnode2]
            full_matrix[subnode1][subnode2] = weight

    full_matrix = np.array(full_matrix)
    return full_matrix


def threshold_matrix(matrix,threshold):
    matrix_size = len(matrix[0])
    new_matrix = create_zero_array(matrix_size)
    new_matrix = np.array(new_matrix)

    for node1 in range(matrix_size):
        for node2 in range(matrix_size):
            weight = matrix[node1][node2]

            if weight >= threshold:
                new_matrix[node1][node2] = weight

    return new_matrix


def threshold_matrix_top_percent(matrix,percent):
    matrix_size = len(matrix[0])

    weight_list = []
    for node1 in range(matrix_size):
        for node2 in range(matrix_size):
            weight = matrix[node1][node2]
            if weight > 0:
                weight_list.append(weight)
    weight_list.sort(reverse=True)


    fraction = float(percent) / 100

    if percent == 0:
        critical_index = 0
    else:
        critical_index = (len(weight_list) * fraction) -1
        critical_index = int(critical_index)
    critical_value = weight_list[critical_index]

    new_matrix = threshold_matrix(matrix,critical_value)
    return new_matrix


def threshold_matrix_window(matrix,minimum,maximum):
    matrix_size = len(matrix[0])

    weight_list = []
    for node1 in range(matrix_size):
        for node2 in range(matrix_size):
            weight = matrix[node1][node2]
            if weight > 0:
                weight_list.append(weight)
    weight_list.sort(reverse=True)

    max_fraction = float(maximum) / 100 #0.1
    min_fraction = float(minimum) / 100 #0.2

    min_critical_index = int((len(weight_list) * min_fraction) - 1)
    min_critical_value = weight_list[min_critical_index]

    if maximum == 0:
        max_critical_index = 0
    else:
        max_critical_index = int((len(weight_list) * max_fraction) - 1)

    max_critical_value = weight_list[max_critical_index]

    new_matrix = create_zero_array(matrix_size)
    new_matrix = np.array(new_matrix)

    for node1 in range(matrix_size):
        for node2 in range(matrix_size):
            weight = matrix[node1][node2]

            if weight >= min_critical_value and weight <= max_critical_value:
                new_matrix[node1][node2] = weight

    return new_matrix




    return new_matrix

def scale_matrix(matrix, new_total):

    original_total = np.sum(matrix)

    if original_total == 0:
        ratio = 0
    else:
        ratio = new_total/original_total

    rows,columns = np.shape(matrix)


    for node_from in range(rows):
        for node_to in range(columns):

            matrix[node_from][node_to] = matrix[node_from][node_to] * ratio

    return matrix


def get_boolean_with_probability(probability):
    number_1 = random.uniform(0,1)

    if probability == 0:
        return False

    if probability == 1:
        return True

    if probability > number_1:
        return True
    else:
        return False


def create_random_binary_matrix(size,density):
    matrix = create_zero_array(size)
    for row in range(size):
        for column in range(size):
            if get_boolean_with_probability(density):
                matrix[row][column] = 1
    matrix = np.array(matrix)
    return matrix


def create_random_binary_list(size,density):
    binary_list = []

    for index in range(size):
        if get_boolean_with_probability(density):
            binary_list.append(1)
        else:
            binary_list.append(0)

    return binary_list



def extract_matrix_region_around_point(matrix,centre_point,radius):

    rows, columns = np.shape(matrix)

    left_border     = int(centre_point[1] - radius)
    right_border    = int((centre_point[1] + radius) + 1)
    top_border      = int(centre_point[0] - radius)
    bottom_border   = int((centre_point[0] + radius) + 1)

    if top_border < 0:
        top_border = 0

    if left_border < 0:
        left_border = 0

    if bottom_border > rows:
        bottom_border = rows

    if right_border > columns:
        right_border = columns

    new_matrix = []

    for row in range(top_border,bottom_border):
       new_matrix.append(matrix[row][left_border:right_border])

    new_matrix = np.array(new_matrix)
    return new_matrix



def convert_matrix_to_tuple(matrix):
    matrix_size, matrix_size_2 = np.shape(matrix)
    new_matrix = create_empty_array(matrix_size)

    for row in range(matrix_size):
        for column in range(matrix_size):
                new_matrix[row][column] = (matrix[row][column], matrix[row][column], matrix[row][column])

    return new_matrix


def get_surrounding_points(point):

    upper_left     = (point[0] - 1, point[1] - 1)
    upper_middle   = (point[0] - 1, point[1] + 0)
    upper_right    = (point[0] - 1, point[1] + 1)
    centre_left    = (point[0] + 0, point[1] - 1)
    centre_right   = (point[0] + 0, point[1] + 1)
    lower_left     = (point[0] + 1, point[1] - 1)
    lower_middle   = (point[0] + 1, point[1] + 0)
    lower_right    = (point[0] + 1, point[1] + 1)

    return [upper_left, upper_middle, upper_right, centre_left, centre_right, lower_left, lower_middle, lower_right]



def get_outline_of_region(matrix, x_points, y_points):
    outline_points = []
    height, width = len(matrix), len(matrix[0])

    #Create A List Of Region Points
    region_points = []
    for index in range(len(x_points)):
        coordinate_pair = (y_points[index],x_points[index])
        region_points.append(coordinate_pair)

    for point in (region_points):
        surronding_points = get_surrounding_points(point)

        for surrounding_point in surronding_points:

            if surrounding_point in region_points:
                pass
            else:
                #if surrounding_point[0] < height and surrounding_point[0] >= 0:
                    #if surrounding_point[1] < width and surrounding_point[1] >= 0:
                 outline_points.append(surrounding_point)


    outline_points = set(outline_points)
    outline_points = list(outline_points)

    return outline_points



def highlight_region(matrix,region,colour):

    for point in region:
        matrix[int(point[0])][int(point[1])] = colour

    return matrix


def circularly_permute(matrix):

    rows, columns = np.shape(matrix)
    permutated_array = []

    for row in range(rows):
        permuted_row = create_empty_list(columns)
        shift = random.randint(0,columns-1)

        for column in range(columns):
            new_column = column + shift
            if new_column >= columns:
                new_column = new_column - columns

            permuted_row[new_column] = matrix[row][column]

        permutated_array.append(permuted_row)
    permutated_array = np.array(permutated_array)
    return permutated_array
