#!/usr/bin/env python3

'''******************************************************************************
 *
 * MAT0122 ALGEBRA LINEAR I
 * Aluno: Gabriel Geraldino de Souza
 * Numero USP: 12543885
 * Tarefa: L02
 * Data: DEZEMBRO, 2021
 * 
 * DECLARO QUE SOU O ÚNICO AUTOR E RESPONSÁVEL POR ESTE PROGRAMA. TODAS AS PARTES 
 * DO PROGRAMA, EXCETO AS QUE SÃO BASEADAS EM MATERIAL FORNECIDO PELO PROFESSOR 
 * FORAM DESENVOLVIDAS POR MIM. DECLARO TAMBÉM QUE SOU RESPONSÁVEL POR TODAS CÓPIAS 
 * DESTE PROGRAMA E QUE NÃO DISTRIBUÍ NEM FACILITEI A DISTRIBUIÇÃO DE CÓPIAS DESTE 
 * PROGRAMA.
 *
 ******************************************************************************'''

import numpy as np
from matplotlib import pyplot as plt

import random
import timeit

def generate_random_matrix(order:int=3) -> list:
	matrix = [[random.random() for i in range(order)] for j in range(order)]
	return matrix

def generate_random_lower_triangular_matrix(order:int=3) -> list:
	matrix = [[0 for i in range(order)] for j in range(order)]
	for i in range(order):
		for j in range(i+1):
			matrix[i][j] = random.random()
	return matrix

def generate_random_upper_triangular_matrix(order:int=3) -> list:
	matrix = [[0 for i in range(order)] for j in range(order)]
	for i in range(order):
		for j in range(order-1, i-1, -1):
			matrix[i][j] = random.random()
	return matrix

def generate_hilbert_matrix(order:int=3) -> list:
	matrix = [[0 for i in range(order)] for j in range(order)]
	for i in range(order):
		for j in range(order):
			matrix[i][j] = 1/(i+j+1)
	return matrix

def inverse_lower(matrix:list) -> list:
	order = len(matrix)
	inverse_matrix = [[0 for i in range(order)] for j in range(order)]
	for i in range(order-1,-1,-1):
		inverse_matrix[i][i] = 1/matrix[i][i]
		for j in range(i+1, order):
			for k in range(i, j):
				inverse_matrix[j][i] += matrix[j][k]*inverse_matrix[k][i]
			inverse_matrix[j][i] /= -matrix[j][j]
	return inverse_matrix

def inverse_upper(matrix:list) -> list:
	order = len(matrix)
	inverse_matrix = [[0 for i in range(order)] for j in range(order)]
	for i in range(order):
		inverse_matrix[i][i] = 1/matrix[i][i]
		for j in range(i+1,order):
			for k in range(i, j):
				inverse_matrix[i][j] += matrix[k][j]*inverse_matrix[i][k]
			inverse_matrix[i][j] /= -matrix[j][j]
	return inverse_matrix

def inverse(matrix:list) -> list:
	order = len(matrix)
	inverse_matrix = [[i for i in row] for row in matrix]
	for i in range(order):
		for j in range(order, 2*order):
			if(j==(i+order)):
				inverse_matrix[i].append(1)
				for _ in range(j+1, 2*order):
					inverse_matrix[i].append(0)
				break
			else:
				inverse_matrix[i].append(0)
	for i in range(order):
		for j in [x for x in range(2*order) if x!=i]:
			inverse_matrix[i][j] /= inverse_matrix[i][i]
		for j in [x for x in range(order) if x!=i]: 
			for k in [x for x in range(2*order) if x!=i]: 
				inverse_matrix[j][k] -= inverse_matrix[j][i] * inverse_matrix[i][k]
	for i in range(order):
		inverse_matrix[i] = inverse_matrix[i][order:]
	return inverse_matrix

def print_matrix(matrix:list) -> None:
	for i in range(len(matrix)):
		for j in range(len(matrix[0])):
			print(round(matrix[i][j],2), end=" ")
		print()

def generate_random_points(qty:int) -> list:
	points = []
	for _ in range(qty):
		x = random.random()*10
		y = random.random()*10
		points.append((x,y))
	return points

def transpose(matrix:list) -> list:
	return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

def multiply(left_matrix:list, right_matrix:list) -> list:
	result = [[0 for i in range(len(right_matrix[0]))] for j in range(len(left_matrix))]
	for i in range(len(left_matrix)):
		for j in range(len(right_matrix[0])):
			for k in range(len(right_matrix)):
				result[i][j] += left_matrix[i][k]*right_matrix[k][j]
	return result

def p(coeffs:list,point:tuple) -> float:
	xi = point[0]
	yi = point[1]
	value = 0
	for i in range(len(coeffs)):
		value += coeffs[i][0]*(xi**i)
	return (value-yi)**2

def deviation(coeffs:list,points:list) -> float:
	total = 0
	for point in points:
		total += p(coeffs,point)
	return total

def polynomial_approximation(points:list,degree:int) -> list:
	A = [[points[j][0]**i for i in range(degree)] for j in range(len(points))]
	b = [[p[1]] for p in points]
	x = multiply(multiply(inverse(multiply(transpose(A), A)), transpose(A)), b)
	return x

def plot_polynomial_approximation(points:list,coeffs:list) -> None:
	title = ""
	for i in range(len(coeffs)):
		title += f"(x^{i})({coeffs[i][0]}) + "
	print(title[:-3])
	plt.title(title[:-3], fontsize=10)
	plt.xlabel('x') 
	plt.ylabel('y') 
	coeffs = np.array(coeffs)
	x = np.linspace(0, 10, 10*len(points))
	y = np.array([np.sum(np.array([coeffs[i]*(j**i) for i in range(len(coeffs))])) for j in x])
	plt.plot(x, y)
	for point in points:
		plt.scatter(np.array([point[0]]), np.array([point[1]]))
	plt.show()
	return

if __name__ == "__main__":
	k = int(input())
	threshold = 0 #useful for a large amount of points; hard to define it when k is small
	points = generate_random_points(k)
	degree = 0
	approximation = polynomial_approximation(points,degree)
	error = deviation(approximation,points)
	while(degree<len(points) and error>threshold):
		degree += 1
		k_approximation = polynomial_approximation(points,degree)
		k_deviation = deviation(k_approximation,points)
		if(k_deviation<error):
			approximation = k_approximation
			error = k_deviation
	print(len(approximation))
	print(error)
	plot_polynomial_approximation(points, approximation)	
