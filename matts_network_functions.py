import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.figure import Figure
from matplotlib import cm
import numpy as np
import random
import math
import community
from scipy import stats


#Functiosn to move to Admin
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


def create_ones_array(size):
    new_array = create_empty_list(size)

    for row in range(size):
        columns = create_list_of(size,float(1))
        new_array[row] = columns

    new_array = np.array(new_array)
    return new_array


def symmetrical_to_upper_triangle(original_matrix):

    number_of_rows,number_of_columns = np.shape(original_matrix)
    new_matrix = create_zero_array(number_of_columns)

    for row in range(number_of_rows):
        for column in range(row, number_of_columns):
            value = original_matrix[row][column]
            new_matrix[row][column] = value

    new_matrix = np.array(new_matrix)
    np.fill_diagonal(new_matrix,0)
    return new_matrix

def create_random_array(size,density,min_value,max_value):
    new_array = create_empty_list(size)

    for row in range(size):
        columns = create_zero_list(size)
        new_array[row] = columns

    new_array = np.array(new_array)

    for row in range(size):
        for column in range(size):
            if get_boolean_with_probability(density) == True:
                value = random.uniform(min_value,max_value)
                new_array[row][column] = value

    return new_array


def create_random_array_int(size,density,min_value,max_value):
    new_array = create_empty_list(size)

    for row in range(size):
        columns = create_zero_list(size)
        new_array[row] = columns

    new_array = np.array(new_array)

    for row in range(size):
        for column in range(size):
            if get_boolean_with_probability(density) == True:
                value = int(random.uniform(min_value,max_value))
                new_array[row][column] = value

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

def create_list_of(size,item):
    new_list = []
    for number in range(size):
        new_list.append(item)

    return new_list


def assign_colour(item,number_of_items,colour_map):
    colour_map = cm.get_cmap(colour_map)
    value = float(item)/number_of_items
    float_tuple = colour_map(value)
    return float_tuple

def create_node_colour_map(node_value_list):

    colour_map = []

    #Get Max Weight
    max_weight = 0
    for node_value_pair in node_value_list:
        weight = node_value_pair[1]
        if weight > max_weight:
            max_weight = weight

    #Assign Colours
    for node_value_pair in node_value_list:
        weight = node_value_pair[1]
        colour = assign_colour(weight,max_weight,"inferno")
        colour_map.append(colour)

    return colour_map



def split_list_into_chunks_of_the_same_thing(input_list):

    unique_values = list(set(input_list))
    unique_values.sort(reverse=True)
    new_list = create_empty_list(len(unique_values))


    for unique_index in range(len(unique_values)):
        unique_item = unique_values[unique_index]

        for list_item in input_list:
            if list_item == unique_item:
                new_list[unique_index].append(list_item)


    return new_list

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



#Module Functions:
def check_module_imported():
    return "Yeppers!"



#Create Graph Functions
def create_graph(connectivity_matrix):
    #graph = nx.from_numpy_matrix(connectivity_matrix,parallel_edges=False)
    graph = nx.from_numpy_array(connectivity_matrix)
    return graph

def create_weighted_graph(connectivity_matrix):
    graph = nx.from_numpy_array(connectivity_matrix, create_using=nx.DiGraph)
    return graph



def get_edge_weights(graph,scaling_factor):
    weight_list = []
    weights = nx.get_edge_attributes(graph,"weight")
    edges = nx.edges(graph)

    for edge in edges:
        weight = weights[edge] * scaling_factor
        weight_list.append(weight)

    return weight_list


def get_edge_weights_thresholded(graph,threshold):
    weight_list = []
    weights = nx.get_edge_attributes(graph, "weight")
    edges = nx.edges(graph)

    for edge in edges:
        weight = weights[edge]

        if weight <= threshold:
            weight_list.append(weight)
        else:
            weight_list.append(0)

    return weight_list



#Clustering Algorithms
def louvain_partition_graph(graph):
    colour_map = []
    partitions = community.best_partition(graph)
    community_list = partitions.values()
    number_of_partitions = max(community_list)

    for neuron in partitions:
        colour_map.append(assign_colour(partitions[neuron],number_of_partitions,"tab10"))#"viridis"

    return colour_map


#Plotting Functions
def plot_graph_basic(graph):
    nx.draw(graph, pos=nx.spring_layout(graph), node_size=40)
    plt.show()

def plot_graph_basic_positions(graph,positions):
    nx.draw(graph, pos=positions, node_size=40)
    plt.show()

def plot_graph_basic_weights(graph):

    edge_weights = get_edge_weights(graph,1)
    colour_map = cm.get_cmap("plasma")
    node_colours = louvain_partition_graph(graph)

    circle_layout = nx.circular_layout(graph)
    #edge_vmin=(min(edge_weights)),edge_vmax=max(edge_weights)
    nx.draw(graph, pos=nx.spring_layout(graph), node_size=40, arrows=False, with_labels=False, edge_color=edge_weights, edge_cmap=colour_map, node_color=node_colours)
    plt.show()

def plot_graph_basic_edge_labels(graph):
    edge_weights = get_edge_weights(graph,0.5)
    #edge_weights = get_edge_weights_thresholded(graph,0.5)
    colour_map = cm.get_cmap("Blues")

    circle_layout = nx.circular_layout(graph)
    #edge_vmin=(min(edge_weights)),edge_vmax=max(edge_weights)
    nx.draw(graph, pos=circle_layout, node_size=40, arrows=False, with_labels=True, edge_color=edge_weights,edge_cmap=colour_map)
    nx.draw_networkx_edge_labels(graph, pos=circle_layout, node_size=40, arrows=False, with_labels=True, edge_color=edge_weights, edge_cmap=colour_map)
    plt.show()


def plot_graph_basic_sizes(graph,sizes):
    edge_weights = get_edge_weights(graph,0.5)
    #edge_weights = get_edge_weights_thresholded(graph,0.5)
    colour_map = cm.get_cmap("plasma")

    circle_layout = nx.circular_layout(graph)
    #edge_vmin=(min(edge_weights)),edge_vmax=max(edge_weights)
    nx.draw(graph, pos=nx.spring_layout(graph), node_size=sizes, arrows=False, with_labels=False, edge_color=edge_weights,edge_cmap=colour_map)
    plt.show()

def plot_graph_basic_sizes_and_positions(graph,sizes,positions,node_colours):
    edge_weights = get_edge_weights(graph,0.5)
    #edge_weights = get_edge_weights_thresholded(graph,0.5)
    colour_map = cm.get_cmap("Purples")

    circle_layout = nx.circular_layout(graph)

    edge_widths = []
    maximum_weight = max(edge_weights)
    for weight in edge_weights:
        edge_widths.append((weight/maximum_weight) * 10)

    #edge_vmin=(min(edge_weights)),edge_vmax=max(edge_weights)
    nx.draw(graph, pos=positions, node_size=40, arrows=False, with_labels=False, edge_color=edge_weights,edge_cmap=colour_map, width=edge_widths, node_color=node_colours)
    plt.show()

def plot_graph_basic_custom_colour_map(graph,positions,node_colours):
    edge_weights = get_edge_weights(graph,0.5)
    colour_map = cm.get_cmap("Purples")

    nx.draw(graph, pos=positions, node_size=200, arrows=False, with_labels=False, edge_color=edge_weights,edge_cmap=colour_map, node_color=node_colours)
    plt.show()



def plot_graph_edge_widths_node_sizes(graph,figure,canvas,node_sizes,node_positions,node_colours):
    edge_weights = get_edge_weights(graph, 0.5)
    colour_map = cm.get_cmap("Purples")

    edge_widths = []
    maximum_weight = max(edge_weights)
    for weight in edge_weights:
        edge_widths.append((weight / maximum_weight) * 10)

    figure.clear()
    axis = figure.add_subplot(111)
    nx.draw(graph, ax=axis, pos=node_positions, node_size=node_sizes, arrows=False, with_labels=False, edge_color=edge_weights,edge_cmap=colour_map, width=edge_widths, node_color=node_colours)
    plt.show()
    canvas.draw()
    canvas.update()


def plot_graph_weights_over_threshold(graph,threshold):

    colour_map = cm.get_cmap("Blues")

    circle_layout = nx.circular_layout(graph)
    threshold_edge_list,edge_weights = get_edges_over_threshold(graph,threshold)


    nx.draw(graph, pos=nx.spring_layout(graph), node_size=40, arrows=False, with_labels=False, edge_color=edge_weights, edge_cmap=colour_map, edge_vmin=(min(edge_weights)),edge_vmax=max(edge_weights),edgelist=threshold_edge_list)
    plt.show()


def plot_graph_custom_colourmap(graph,figure,canvas,scaled_node_positions,custom_colour_map,threshold,node_sizes=40):
    figure.clear()
    axis = figure.add_subplot(111)
    colour_map = cm.get_cmap("Blues")
    threshold_edge_list, edge_weights = get_edges_over_threshold(graph, threshold)
    node_colours = custom_colour_map

    nx.draw(graph, pos=scaled_node_positions, ax=axis, node_size=node_sizes, arrows=False, with_labels=False, edge_color=edge_weights,edge_cmap=colour_map,node_color = node_colours, edgelist=threshold_edge_list, edge_vmax=max(edge_weights), edge_vmin=0)
    plt.show()
    canvas.draw()
    canvas.update()

def plot_graph_custom_node_and_edge_colourmaps(graph,figure,canvas,scaled_node_positions,node_colours, edge_colours):
    figure.clear()
    axis = figure.add_subplot(111)
    edge_weights = get_edge_weights(graph,1)
    nx.draw(graph, pos=scaled_node_positions, ax=axis, node_size=40, arrows=False, with_labels=False, edge_color=edge_colours,node_color = node_colours, edge_vmax=max(edge_weights), edge_vmin=1)
    plt.show()
    canvas.draw()
    canvas.update()

def plot_graph_custom_colourmap_sizes(graph,figure,canvas,scaled_node_positions,custom_colour_map,threshold,node_sizes):
    figure.clear()
    axis = figure.add_subplot(111)
    colour_map = cm.get_cmap("Blues")
    threshold_edge_list, edge_weights = get_edges_over_threshold(graph, threshold)
    node_colours = custom_colour_map

    nx.draw(graph, pos=scaled_node_positions, ax=axis, node_size=node_sizes, arrows=False, with_labels=False,
            edge_color=edge_weights, edge_cmap=colour_map, node_color=node_colours, edgelist=threshold_edge_list)
    plt.show()
    canvas.draw()
    canvas.update()


def plot_graph_thresholed_weights(figure,canvas,graph,scaled_node_positions,threshold):
    figure.clear()
    axis = figure.add_subplot(111)
    colour_map = cm.get_cmap("Blues")
    threshold_edge_list, edge_weights = get_edges_over_threshold(graph, threshold)
    node_colours = louvain_partition_graph(graph)
    node_sizes = create_list_of(len(scaled_node_positions),40)#get_centrality_node_sizes(graph,750)
    nx.draw(graph, pos=scaled_node_positions, ax=axis, node_size=40, arrows=False, with_labels=False, edge_color=edge_weights,edge_cmap=colour_map, edgelist=threshold_edge_list,node_color = node_colours)
    plt.show()
    canvas.draw()
    canvas.update()


def plot_graph_thresholed_weights_edge_colour_map(figure,canvas,graph,scaled_node_positions,threshold,edge_colour_map):
    figure.clear()
    axis = figure.add_subplot(111)
    colour_map = cm.get_cmap(edge_colour_map)
    threshold_edge_list, edge_weights = get_edges_over_threshold(graph, threshold)
    node_colours = louvain_partition_graph(graph)
    node_sizes = get_centrality_node_sizes(graph,750)
    nx.draw(graph, pos=scaled_node_positions, ax=axis, node_size=40, arrows=False, with_labels=False, edge_color=edge_weights,edge_cmap=colour_map, edgelist=threshold_edge_list,node_color = node_colours)
    plt.show()
    canvas.draw()
    canvas.update()

def plot_graph_node_sizes(figure,canvas,graph,edge_weights,node_sizes,scaled_node_positions):
    spring_layout = nx.spring_layout(graph)
    circle_layout = nx.circular_layout(graph)
    kamada_layout = nx.kamada_kawai_layout(graph)
    spectral_layout = nx.spectral_layout(graph)

    figure.clear()
    axis = figure.add_subplot(111)

    colour_map = cm.get_cmap("Blues")
    node_colours = louvain_partition_graph(graph)

    nx.draw(graph,pos=scaled_node_positions, ax=axis, node_size=node_sizes,node_color = node_colours, edge_color=edge_weights, edge_vmin=0, edge_vmax=1, edge_cmap=colour_map, with_labels=False, font_color="white", font_size="10")

    canvas.draw()
    canvas.update()


def default_plot_graph(figure,canvas,graph,positions=None,node_sizes=40,colour_map="Blues",node_colours="b",edge_weights=None,edge_scaling_factor=1):
    figure.clear()
    axis = figure.add_subplot(111)
    map = cm.get_cmap(colour_map)

    if positions == None:
        positions = nx.spring_layout(graph)

    if edge_weights == None:
        edge_weights = get_edge_weights(graph,edge_scaling_factor)

    nx.draw(graph,
            ax=axis,
            pos=positions,
            node_size=node_sizes,
            node_color=node_colours,
            edge_color=edge_weights,
            edge_vmin=0,
            edge_vmax=1,
            edge_cmap=map)

    canvas.draw()
    canvas.update()


def get_centrality_node_sizes(graph,scaling_factor):
    eigenvector_centrailities = nx.algorithms.eigenvector_centrality(graph)
    number_of_nodes = len(eigenvector_centrailities)
    centraility_list = create_empty_list(number_of_nodes)

    for node in range(number_of_nodes):
        centraility_list[node] = (eigenvector_centrailities[node] * eigenvector_centrailities[node]) * scaling_factor * 5

    #print "centraility list", centraility_list
    return centraility_list


def expanded_subnetwork_analysis(graphs, original_subnetwork):
    expanded_subnetworks = []
    expanded_node_dict = {}

    for graph in graphs:
        expanded_subnetwork = get_expanded_subnetwork(graph,original_subnetwork)
        expanded_subnetworks.append(expanded_subnetwork)

    number_of_subnetworks = len(expanded_subnetworks)
    for subnetwork in expanded_subnetworks:
        for node in subnetwork:
            if node in expanded_node_dict:
                number_of_appearences = expanded_node_dict[node]
                expanded_node_dict[node] = number_of_appearences + 1
            else:
                expanded_node_dict[node] = 1


    for node in expanded_node_dict.keys():
        expanded_node_dict[node] = (float(expanded_node_dict[node])/number_of_subnetworks) * 100

    return expanded_node_dict




    #return expanded_subnetworks



def get_expanded_subnetwork(graph,original_subnetwork):
    additional_nodes = []
    edge_dict = nx.get_edge_attributes(graph,"weight")

    for node in original_subnetwork:
        strongest_neighbour = None
        strongest_weight = 0

        neighbours = list(nx.neighbors(graph,node))

        for neighbour in neighbours:

            if neighbour in original_subnetwork:
                pass
            else:
                if (node,neighbour) in edge_dict:
                    edge_weight = edge_dict[(node,neighbour)]

                elif (neighbour,node) in edge_dict:
                    edge_weight = edge_dict[(neighbour,node)]

                if edge_weight > strongest_weight:
                        strongest_weight = edge_weight
                        strongest_neighbour = neighbour

        additional_nodes.append(strongest_neighbour)

    additional_nodes = set(additional_nodes)
    additional_nodes = list(additional_nodes)
    additional_nodes.sort()

    return additional_nodes
    #print additional_nodes








def convert_matrix_to_binary(matrix):
    rows,columns = np.shape(matrix)
    new_matrix = create_zero_array(rows)

    for row in range(rows):
        for column in range(columns):
            if matrix[row][column] > 0:
                new_matrix[row][column] = 1

    return new_matrix




#Rewire Functions
def randomly_rewire(graph,probability):
    edges = nx.edges(graph)
    number_of_nodes = nx.number_of_nodes(graph)

    for edge in edges:

        #Rewire with certain probability
        if get_boolean_with_probability(probability) == True:

            node_1 = edge[0]
            node_2 = edge[1]

            #Create a list of possible new nodes, avoiding existing connections
            existing_neighbours = nx.neighbors(graph,node_1)
            list_of_new_potential_neighbours = range(0,number_of_nodes)
            for neighbour in existing_neighbours:
                list_of_new_potential_neighbours.remove(neighbour)

            #Pick the new random node
            new_node_index = random.randint(0,(len(list_of_new_potential_neighbours)-1))
            new_node = list_of_new_potential_neighbours[new_node_index]

            #Remove old connection and add the new one
            graph.remove_edge(node_1,node_2)
            graph.add_edge(node_1, new_node)

    return graph



def create_list_of_integers(size):
    new_list = []
    for x in range(size):
        list.append(x)

    return new_list








def get_upper_triangle_part_of_matrix(matrix):

    list_of_items = []
    number_of_rows, number_of_columns = np.shape(matrix)

    for row in range(number_of_rows-1):
        for column in range(row+1, number_of_columns):
            item = matrix[row][column]
            list_of_items.append(item)

    return list_of_items


def randomise_graph(graph):

    number_of_nodes = nx.number_of_nodes(graph)
    random_matrix = create_zero_array(number_of_nodes)
    connection_matrix = nx.to_numpy_array(graph)

    # Create a list of randomised edges
    randomised_edges = get_upper_triangle_part_of_matrix(connection_matrix)
    random.shuffle(randomised_edges)

    # Go Through Upper Triangular Matrix and Fill in Edges
    edge = 0
    for row in range(number_of_nodes-1):
        for column in range(row+1, number_of_nodes):
            random_matrix[row][column] = randomised_edges[edge]
            edge += 1

    # Turn Matrix into Graph
    random_graph = create_graph(random_matrix)
    return random_graph




def randomly_rewire_undirected(graph,probability):

    number_of_rewirings = 0
    number_of_nodes = nx.number_of_nodes(graph)
    weight_dict = nx.get_edge_attributes(graph, "weight")

    existing_edges = nx.edges(graph)
    node_list = nx.nodes(graph)

    #Get List Of All Potential Edges (undirected):
    all_potential_edges = []
    for node_1 in node_list:
        for node_2 in node_list:
            if node_1 == node_2:
                pass
            else:
                all_potential_edges.append((node_1,node_2))

    #Get Pool of Re-wired Edges
    potential_pool = all_potential_edges
    for edge in existing_edges:

        if edge in potential_pool: #Should I have Had to Add This In?
            potential_pool.remove(edge)

        if (edge[1],edge[0]) in potential_pool:
            potential_pool.remove((edge[1],edge[0]))

    #Go Through and Poentially Swap Each One
    new_graph_edge_dict = {}

    for edge in existing_edges:

        #If Edge Not Swapped Just Add To The New Edge List
        if get_boolean_with_probability(probability) == False or len(potential_pool) == 0:
            new_graph_edge_dict[edge] = weight_dict[edge]

        else:
            number_of_rewirings += 1
            #Select an Edge to swap with
            new_edge = potential_pool[random.randint(0,len(potential_pool)-1)]

            #Add New Edge To New Graph
            new_graph_edge_dict[new_edge] = weight_dict[edge]

            #Remove New Edge From Potential Pool
            potential_pool.remove(new_edge)

            #Add Previous Edge To Potential Pool
            potential_pool.append(edge)


    #Turn New Edge Dict Into New Graph
    new_graph = nx.Graph()
    for node in node_list:
        new_graph.add_node(node)

    for edge in new_graph_edge_dict:
        new_graph.add_edge(edge[0], edge[1], weight=new_graph_edge_dict[edge])

    return new_graph



#Lattice Graph Functions
def create_directed_unweighted_lattice(number_of_nodes):
    connection_matrix = create_zero_array(number_of_nodes)

    for node_from in range(number_of_nodes):
        for neighbour in range(-2,3):

            if neighbour == 0:
                pass
            else:

                node_to = node_from + neighbour

                if node_to < 0:
                    node_to = number_of_nodes + node_to

                if node_to >= (number_of_nodes):
                    node_to = node_to - (number_of_nodes + 1)

                connection_matrix[node_from][node_to] = 1

    np.fill_diagonal(connection_matrix, 0)
    graph = create_graph(connection_matrix)
    return graph

def create_directed_weighted_lattice(number_of_nodes):
    connection_matrix = create_zero_array(number_of_nodes)
    spacing = 2  # + 1
    start_value = 0 - spacing
    end_value = spacing + 1

    for node_from in range(number_of_nodes):
        for neighbour in range(start_value, end_value):

            if neighbour == 0:
                pass
            else:
                node_to = node_from + neighbour

                if node_to == node_from:
                    pass

                else:

                    if node_to < 0:
                        node_to = (number_of_nodes+node_from) + neighbour

                    elif node_to >= (number_of_nodes):
                        node_to = node_to - (number_of_nodes)


                    maximum_distance = (float(number_of_nodes) / 2) + 1


                    current_distance = abs(neighbour)
                    weight = maximum_distance - current_distance

                    connection_matrix[node_from][node_to] = weight

    np.fill_diagonal(connection_matrix, 0)

    graph = create_weighted_graph(connection_matrix)
    return graph

def create_undirected_weighted_lattice(number_of_nodes):
    connection_matrix = create_zero_array(number_of_nodes)
    spacing = 2  # + 1
    start_value = 0 - spacing
    end_value = spacing + 1

    for node_from in range(number_of_nodes):
        for neighbour in range(start_value, end_value):

            if neighbour == 0:
                pass
            else:
                node_to = node_from + neighbour

                if node_to == node_from:
                    pass

                else:

                    if node_to < 0:
                        node_to = (number_of_nodes+node_from) + neighbour

                    elif node_to >= (number_of_nodes):
                        node_to = node_to - (number_of_nodes)


                    maximum_distance = (float(number_of_nodes) / 2) + 1


                    current_distance = abs(neighbour)
                    weight = maximum_distance - current_distance

                    if connection_matrix[node_to][node_from] == 0:
                        connection_matrix[node_from][node_to] = weight

    np.fill_diagonal(connection_matrix, 0)

    graph = create_weighted_graph(connection_matrix)
    return graph





def create_directed_equivalent_lattice(graph):
    connection_matrix = nx.to_numpy_array(graph)
    number_of_nodes, number_of_nodes_2 = np.shape(connection_matrix)
    edge_weights = get_edge_weights(graph, 1)
    weight_groups = split_list_into_chunks_of_the_same_thing(edge_weights)
    lattice_matrix = create_zero_array(number_of_nodes)

    distance_index = 1
    points_to_fill = get_matrix_positions_at_distance(lattice_matrix, distance_index)

    for group in weight_groups:
        for weight in group:

            # Fill that point and remove it from pool
            point_index = random.randint(0, len(points_to_fill) - 1)
            point_to_fill = points_to_fill[point_index]
            lattice_matrix[point_to_fill[0]][point_to_fill[1]] = weight

            # Remove the reverse point to make sure it dosent take a different weight
            points_to_fill.remove(point_to_fill)

            # If the list is empty, get the next distance
            if len(points_to_fill) == 0:
                distance_index += 1
                new_potential_points = get_matrix_positions_at_distance(lattice_matrix, distance_index)

                # Check They are all empty and their reverse are not filled
                for new_potential_point in new_potential_points:
                    if lattice_matrix[new_potential_point[0]][new_potential_point[1]] == 0:
                        points_to_fill.append(new_potential_point)

    lattice_matrix = np.array(lattice_matrix)
    lattice_graph = create_weighted_graph(lattice_matrix)
    return lattice_graph

def create_undirected_equivalent_lattice(graph):

    #graph = nx.to_undirected(graph)

    connection_matrix = nx.to_numpy_array(graph)
    number_of_nodes, number_of_nodes_2 = np.shape(connection_matrix)
    edge_weights = get_edge_weights(graph, 1)
    weight_groups = split_list_into_chunks_of_the_same_thing(edge_weights)
    lattice_matrix = create_zero_array(number_of_nodes)

    distance_index = 1
    points_to_fill = get_undirected_matrix_positions_at_distance(lattice_matrix, distance_index)

    for group in weight_groups:
        for weight in group:

            # Fill that point and remove it from pool
            point_index = random.randint(0, len(points_to_fill) - 1)


            point_to_fill = points_to_fill[point_index]
            lattice_matrix[point_to_fill[0]][point_to_fill[1]] = weight


            # Remove the reverse point to make sure it dosent take a different weight
            points_to_fill.remove(point_to_fill)

            # If the list is empty, get the next distance
            if len(points_to_fill) == 0:
                distance_index += 1
                new_potential_points = get_undirected_matrix_positions_at_distance(lattice_matrix, distance_index)

                # Check They are all empty and their reverse are not filled
                for new_potential_point in new_potential_points:
                    if lattice_matrix[new_potential_point[0]][new_potential_point[1]] == 0:
                        points_to_fill.append(new_potential_point)

    lattice_matrix = np.array(lattice_matrix)
    lattice_graph = create_weighted_graph(lattice_matrix)
    return lattice_graph



#Get Matrix Positions
def get_matrix_positions_at_distance(matrix,distance):
    matrix_height,matrix_width = np.shape(matrix)

    list_of_positions = []

    for row in range(matrix_height):
            node_in_front = row + distance
            node_behind = row - distance

            if node_in_front >= (matrix_width):
                node_in_front = node_in_front - (matrix_width)

            if node_behind < 0:
                node_behind = matrix_width + node_behind

            if node_in_front == node_behind:
                list_of_positions.append([row,node_in_front])

            else:
                list_of_positions.append((row, int(node_behind)))
                list_of_positions.append((row, int(node_in_front)))

    return list_of_positions

def get_undirected_matrix_positions_at_distance(matrix,distance):
    matrix_height,matrix_width = np.shape(matrix)

    list_of_positions = []

    for row in range(matrix_height):
            node_in_front = row + distance
            node_behind = row - distance

            if node_in_front >= (matrix_width):
                node_in_front = node_in_front - (matrix_width)

            if node_behind < 0:
                node_behind = matrix_width + node_behind

            if node_in_front == node_behind:
                list_of_positions.append([row,node_in_front])

            else:
                list_of_positions.append((row, int(node_behind)))
                list_of_positions.append((row, int(node_in_front)))

    for position in list_of_positions:
        reverse_position = (position[1],position[0])
        if reverse_position in list_of_positions:
            list_of_positions.remove((reverse_position))

    return list_of_positions


def get_edges_over_threshold(graph,threshold):

    edge_list = []
    weight_list = []

    edges = nx.get_edge_attributes(graph,"weight")

    for edge in edges:
        if edges[edge] >= threshold:
            edge_list.append(edge)
            weight_list.append(edges[edge])

    return edge_list, weight_list


def invert_weights(graph):
    connection_matrix = nx.to_numpy_array(graph)
    #connection_matrix = symmetrical_to_upper_triangle(connection_matrix)


    number_of_nodes = nx.number_of_nodes(graph)

    for row in range(number_of_nodes):
        for column in range(number_of_nodes):
            value = connection_matrix[row][column]
            if value == 0:
                pass
            else:
                value = float(1/value)
                connection_matrix[row][column] = value

    graph = create_graph(connection_matrix)
    return graph


def invert_weights_and_scale(graph):

    connection_matrix = nx.to_numpy_array(graph)
    connection_matrix = symmetrical_to_upper_triangle(connection_matrix)

    number_of_nodes = nx.number_of_nodes(graph)

    for row in range(number_of_nodes):
        for column in range(number_of_nodes):
            value = connection_matrix[row][column]
            if value == 0:
                pass
            else:
                value = float(1) / float(value)
                connection_matrix[row][column] = value

    graph = create_graph(connection_matrix)
    return graph


#Graph Analysis Functions
def get_small_world_propensity(graph):

    equivalent_lattice = create_undirected_equivalent_lattice(graph)
    equivalent_random = randomise_graph(graph)

    observed_clustering = nx.average_clustering(graph,weight="weight")
    lattice_clustering  = nx.average_clustering(equivalent_lattice,weight="weight")
    random_clustering   = nx.average_clustering(equivalent_random,weight="weight")

    #plot_graph_basic_weights(graph)
    #plot_graph_basic_weights(equivalent_lattice)
    #plot_graph_basic_weights(equivalent_random)

    observed_length = get_average_shortest_path_length(graph)
    lattice_length  = get_average_shortest_path_length(equivalent_lattice)
    random_length   = get_average_shortest_path_length(equivalent_random)

    try:
        delta_c = (lattice_clustering - observed_clustering) / (lattice_clustering - random_clustering)
        delta_l = (observed_length - random_length) / (lattice_length - random_length)
        small_world_propensity = math.sqrt(((delta_c*delta_c) + (delta_l*delta_l)) / 2)
    except:
        small_world_propensity = 0
    """
    print "observed_clustering" , observed_clustering
    print "lattice_clustering"  , lattice_clustering
    print "random_clustering"   , random_clustering

    print "observed_length"     , observed_length
    print "lattice_length"      , lattice_length
    print "random_length"       , random_length

    print "delta c", delta_c
    print "delta l", delta_l

    print "delta c sqaured", delta_c * delta_c
    print "delta l sqaured", delta_l * delta_l
    print "sum ", (delta_c*delta_c) + (delta_l*delta_l)
    print "halved", ((delta_c*delta_c) + (delta_l*delta_l)) / 2
    print "small world propensity", small_world_propensity
    """

    return small_world_propensity


def get_small_world_index(graph):

    equivalent_random = randomise_graph(graph)

    observed_clustering = nx.average_clustering(graph, weight="weight")
    random_clustering = nx.average_clustering(equivalent_random, weight="weight")

    observed_length = get_average_shortest_path_length(graph)
    random_length = get_average_shortest_path_length(equivalent_random)

    try:
        small_world_index = (observed_clustering/random_clustering) / (observed_length/random_length)
    except:
        small_world_index = 0

    return small_world_index

def get_average_degree(graph):
    degree_list = []
    node_degree_pairs = list(nx.degree(graph, weight="weight"))

    for pair in node_degree_pairs:
        degree_list.append(pair[1])

    average_degree = np.mean(degree_list)
    return average_degree

def get_degree_variability(graph):
    degree_list = []
    node_degree_pairs = list(nx.degree(graph,weight="weight"))

    for pair in node_degree_pairs:
        degree_list.append(pair[1])

    degree_variability = np.var(degree_list)

    return degree_variability


def get_average_shortest_path_length(graph):

    inverse_graph = invert_weights(graph)
    inverse_graph = nx.to_undirected(inverse_graph)

    try:
        shortest_path_length = nx.average_shortest_path_length(inverse_graph, weight="weight")
        shortest_path_length = 1 / float(shortest_path_length)
    except:
        shortest_path_length = 999

    return shortest_path_length


def get_average_graph_weight(graph):

    number_of_edges = nx.number_of_edges(graph)
    total_weight = 0

    for from_node in graph.nodes():
        for to_node in graph.nodes():
            weight_dict = graph.get_edge_data(from_node, to_node)
            if weight_dict == None:
                pass
            else:
                weight = weight_dict["weight"]
                total_weight += weight

    average_weight = total_weight / (number_of_edges*2)
    return average_weight


def get_louvain_clusters(graph):

    #Parition Graph Into Clusters Using Louvain Algorithm
    partitions = community.best_partition(graph)

    #Create Lists of the Nodes in Each Component
    community_list = partitions.values()
    number_of_clusters = max(community_list) + 1
    subgraph_node_lists = create_empty_list(number_of_clusters)
    for node in range(len(community_list)):
        node_community = community_list[node]
        subgraph_node_lists[node_community].append(node)

    return subgraph_node_lists


def perform_cluster_analysis(graph):

    #Parition Graph Into Clusters Using Louvain Algorithm
    partitions = community.best_partition(graph)

    #Create Lists of the Nodes in Each Component
    community_list = partitions.values()
    number_of_clusters = max(community_list) + 1
    subgraph_node_lists = create_empty_list(number_of_clusters)
    for node in range(len(community_list)):
        node_community = community_list[node]
        subgraph_node_lists[node_community].append(node)

    #Get the within Module Strengths
    within_cluster_average_strength = []
    subgraphs = []
    edges = []
    for subgraph_node_list in subgraph_node_lists:
        subgraph = nx.subgraph(graph, subgraph_node_list)
        subgraphs.append(subgraph)
        average_weight = get_average_graph_weight(subgraph)
        within_cluster_average_strength.append(average_weight)

    #Get The Between Module Strengths
    connectivity_matrix = create_zero_array(number_of_clusters)
    for cluster_from in range(number_of_clusters):
        for cluster_to in range(number_of_clusters):
            if cluster_from == cluster_to:
                pass
            else:
                total_weight = 0
                for node_from in subgraphs[cluster_from]:
                    for node_to in subgraphs[cluster_to]:
                        weight_dict = graph.get_edge_data(node_from, node_to)
                        if weight_dict == None:
                            pass
                        else:
                            weight = weight_dict["weight"]
                            total_weight += weight

                total_weight = total_weight / (len(subgraphs[cluster_from]) + len(subgraphs[cluster_to]))
                connectivity_matrix[cluster_from][cluster_to] = total_weight


    return number_of_clusters, subgraph_node_lists, within_cluster_average_strength, connectivity_matrix


def create_shuffled_conditions(list_1,list_2):

    shuffle_condition_1 = create_empty_list(len(list_1))
    shuffle_condition_2 = create_empty_list(len(list_2))
    combined_list = list_1 + list_2

    for item_space in range(len(shuffle_condition_1)):
        item_index = random.randint(0,(len(combined_list)-1))
        item = combined_list[item_index]
        shuffle_condition_1[item_space] = item
        del combined_list[item_index]

    shuffle_condition_2 = combined_list

    return shuffle_condition_1, shuffle_condition_2



def detect_edge_differences_in_graph_sets(matrix_set_1,matrix_set_2):

    t_stat_matrix, p_value_matrix = create_edge_difference_matrix(matrix_set_1,matrix_set_2)
    permutation_threshold_matrix = get_edge_difference_permutation_statistics(matrix_set_1,matrix_set_2,200)

    #Return Only Signficant T Stats
    number_of_nodes,number_of_nodes_2 = np.shape(matrix_set_1[0])
    signficant_t_stats = create_zero_array(number_of_nodes)

    for from_node in range(number_of_nodes):
        for to_node in range(number_of_nodes):

            p_value = p_value_matrix[from_node][to_node]

            if p_value < 0.05:
                threshold = permutation_threshold_matrix[from_node][to_node]
                t_stat = t_stat_matrix[from_node][to_node]
                abs_t_stat = abs(t_stat)
                if abs_t_stat > threshold:
                    signficant_t_stats[from_node][to_node] = t_stat

    return signficant_t_stats


def get_edge_difference_permutation_statistics(matrix_set_1,matrix_set_2,number_of_permutations):

    number_of_nodes,number_of_nodes_2 = np.shape(matrix_set_1[0])

    permutation_list_matrix = create_empty_array(number_of_nodes)
    permutation_threshold_matrix = create_empty_array(number_of_nodes)

    for permutation in range(number_of_permutations):
        print "permutation: ",permutation

        # First Create Shuffled Sample
        shuffle_condition_1, shuffle_condition_2 = create_shuffled_conditions(matrix_set_1,matrix_set_2)

        #Get T and P Matricies for Shuffled Condition
        t_stat_matrix, p_value_matrix = create_edge_difference_matrix(shuffle_condition_1,shuffle_condition_2)

        #Add Each Value To Permutation List Matrix
        for from_node in range(number_of_nodes):
            for to_node in range(number_of_nodes):
                stat = t_stat_matrix[from_node][to_node]
                permutation_list_matrix[from_node][to_node].append(abs(stat))


    #Once Through All The Permutations, Add The Threshold To The Threshold Matrix
    for from_node in range(number_of_nodes):
        for to_node in range(number_of_nodes):
            permutation_list = permutation_list_matrix[from_node][to_node]
            permutation_list.sort()
            threshold_index = int((len(permutation_list) * 0.95) - 1)
            threshold_value = permutation_list[threshold_index]
            permutation_threshold_matrix[from_node][to_node] = threshold_value

    return permutation_threshold_matrix


def create_edge_difference_matrix(connection_matrix_list_1,connection_matrix_list_2):

    number_of_nodes, number_of_nodes_2 = np.shape(connection_matrix_list_1[0])

    p_value_matrix = create_ones_array(number_of_nodes)
    t_stat_matrix = create_zero_array(number_of_nodes)

    for node in range(number_of_nodes):
        t_stats, p_values = detect_edge_differences(node,connection_matrix_list_1,connection_matrix_list_2,number_of_nodes)

        for stat_index in range(len(t_stats)):
            stat = t_stats[stat_index]
            if math.isnan(stat):
                pass
            else:
                t_stat_matrix[node][stat_index] = stat

        for p_value_index in range(len(p_values)):
            p_value = p_values[p_value_index]
            if math.isnan(p_value):
                pass
            else:
                p_value_matrix[node][p_value_index] = p_value

    return t_stat_matrix, p_value_matrix


def detect_edge_differences(source_node,condition_1_connection_matricies,condition_2_connection_matricies,number_of_nodes):

    #Get Edges For Condition 1
    condition_1_weight_list = create_empty_list(number_of_nodes)
    for matrix in condition_1_connection_matricies:
        for node in range(number_of_nodes):
            condition_1_weight_list[node].append(matrix[source_node][node])

    #Get Edges For Condition 2
    condition_2_weight_list = create_empty_list(number_of_nodes)
    for matrix in condition_2_connection_matricies:
        for node in range(number_of_nodes):
            condition_2_weight_list[node].append(matrix[source_node][node])

    #Perform T-Tests
    t_stats = []
    p_values = []

    for node in range(number_of_nodes):
        if node == source_node:
            t_stats.append(0)
            p_values.append(1)
        else:
            t_stat, p_value = stats.ttest_ind(condition_1_weight_list[node],condition_2_weight_list[node])
            t_stats.append(t_stat)
            p_values.append(p_value)

    return t_stats, p_values


def average_eigenvector_centrailtiy(graph):
    average_centraility = 0

    try:
        centrality = nx.eigenvector_centrality(graph,weight="weight",max_iter=10000)
        average_centraility = np.mean(centrality.values())
        average_centraility = np.around(average_centraility,6)
    except:
        pass

    return average_centraility

def false_discovery_rate(p_values,alpha):

    p_values.sort()
    number_of_tests = len(p_values)
    q_values = []
    decision_list = []

    #Get List of Q Values
    for index in range(number_of_tests):
        q_value = alpha * (float((index+1))/number_of_tests)
        q_values.append(q_value)

    #Threshold Each P Value by Q Value
    for index in range(number_of_tests):
        p_value = p_values[index]
        q_value = q_values[index]

        if p_value < q_value:
            decision_list.append(1)
        else:
            decision_list.append(0)

    #Find Highest P Value That Passes Threshold
    highest_passing_p = -1
    for index in range(len(decision_list)):
        if decision_list[index] == 1:
            highest_passing_p = index

    #Set all P values below highest passing p to True
    for index in range(highest_passing_p+1):
        decision_list[index] = 1

    print p_values
    print q_values
    return decision_list






"""
random_matrix = create_random_array(4,0.5,1,10)
random_matrix = symmetrical_to_upper_triangle(random_matrix)
random_graph = create_graph(random_matrix)

plot_graph_basic_weights(random_graph)

inverted_graph = invert_weights(random_graph)
plot_graph_basic_weights(inverted_graph)
plot_graph_basic_weights(random_graph)

#matrix = create_random_array(5,0.9,0,0.1)
#matrix = np.array(matrix)
#np.fill_diagonal(matrix,0)



matrix = [
            [0,0,3,2,0],
            [0,0,4,6,5],
            [3,4,0,3,0],
            [2,6,3,0,1],
            [0,5,0,1,0]
        ]

matrix = np.array(matrix)


matrix = create_random_array_int(10,0.4,0,11)
matrix = symmetrical_to_upper_triangle(matrix)
np.fill_diagonal(matrix,0)
graph = create_graph(matrix)
number_of_clusters, subgraph_node_lists, within_cluster_average_strength = perform_cluster_analysis(graph)

print "number of clusters", number_of_clusters
print "node lsits", subgraph_node_lists
print "within average strengths", within_cluster_average_strength

plot_graph_basic_edge_labels(graph)"""

"""
condition_1_list = [
[[0, 1.1, 1.2], [1.0, 0, 5.0], [1.1, 5.0, 0]],
[[0, 1.0, 1.2], [1.0, 0, 6.0], [1.1, 6.0, 0]],
[[0, 1.1, 1.1], [1.1, 0, 7.0], [1.2, 7.0, 0]],
[[0, 1.1, 1.1], [1.2, 0, 8.0], [1.0, 8.0, 0]],
[[0, 1.0, 1.2], [1.3, 0, 9.0], [1.3, 9.0, 0]],
]

condition_2_list = [
[[0, 1.1, 1.1], [1.1, 0, 0.5], [1.2, 0.5, 0]],
[[0, 1.2, 1.2], [1.2, 0, 0.4], [1.0, 0.4, 0]],
[[0, 1.0, 1.0], [1.1, 0, 0.3], [1.1, 0.3, 0]],
[[0, 1.0, 1.1], [1.2, 0, 0.2], [1.3, 0.2, 0]],
[[0, 1.3, 1.2], [1.0, 0, 0.1], [1.2, 0.1, 0]],
]

t_stat_matrix, p_value_matrix = detect_edge_differences_in_graph_sets(condition_1_list,condition_2_list)
print "P Value Matrix", p_value_matrix

# Turn P Value Matrix Into Connectome
p_weight_matrix = create_empty_array(3)
for from_node in range(3):
    for to_node in range(3):
        p_value = p_value_matrix[from_node][to_node]

        if p_value < 0.05:
            p_weight_matrix[from_node][to_node] = 10  # - p_value
            print from_node, to_node
        else:
            p_weight_matrix[from_node][to_node] = 0

p_weight_matrix = np.array(p_weight_matrix)
print p_weight_matrix
edge_difference_graph = create_graph(p_weight_matrix)
plot_graph_basic_weights(edge_difference_graph)
"""