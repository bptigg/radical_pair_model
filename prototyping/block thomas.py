import numpy as np

matrix = [[1,2,3,1,0,0,0,0,0],
          [4,0,6,0,1,0,0,0,0],
          [7,8,9,0,0,1,0,0,0],
          [1,0,0,9,8,7,1,0,0],
          [0,1,0,6,0,4,0,1,0],
          [0,0,1,3,2,1,0,0,1],
          [0,0,0,1,0,0,1,2,3],
          [0,0,0,0,1,0,4,0,6],
          [0,0,0,0,0,1,7,8,9]]

results = [3,6,3,2,7,2,9,1,9]

block_size = 3
matrix_size = 9
step = int(matrix_size / block_size)

blocks = []
result_blocks = []
for i in range(0,step):
    result_blocks.append([])
    for e in range(0,step):
        result_blocks[i].append(results[(i*block_size + e)])
        if(e == i - 1 or e == i or e == i + 1):
            block = []
            size = block_size**2
            row = -1
            offset_y = i * block_size
            offset_x = e * block_size
            for j in range(0,size):
                if(j % block_size == 0):
                    block.append([])
                    row = row + 1
                block[row].append(matrix[offset_y + row][offset_x + (j%block_size)])
            blocks.append(block)

print(blocks)

row_blocks = []
new_results = []

for i in range(0,step):
    row_blocks.append([])
    #new_results.append()
    e = (3 * i) - 1
    if(i == 0):
        a = np.array(blocks[0])
        ainv = np.linalg.inv(a)
        row_blocks[i].append(np.identity(block_size))
        row_blocks[i].append(ainv @ blocks[1])
        new_results.append(ainv @ result_blocks[i])
    elif(i == step - 1):
        a = np.array(blocks[e+1]- np.array(blocks[e]) @ row_blocks[i-1][1])
        ainv = np.linalg.inv(a)
        row_blocks[i].append(np.identity(block_size))
        new_results.append(ainv @ (result_blocks[i] - np.array(blocks[e]) @ new_results[i-1]))
    else:
        a = np.array(blocks[e+1]- np.array(blocks[e]) @ row_blocks[i-1][1])
        ainv = np.linalg.inv(a)
        row_blocks[i].append(np.identity(block_size))
        row_blocks[i].append(ainv @ blocks[e+2])
        new_results.append(ainv @ (result_blocks[i] - np.array(blocks[e]) @ new_results[i-1]))

print(row_blocks)
print(new_results)

x_vec = []
for i in range(0,matrix_size):
    x_vec.append(0)

for e in range(0,len(new_results[step-1])):
        x_vec[matrix_size - len(new_results[step-1]) + e] = new_results[step-1][e]

num = 2
for i in range(step-2, -1, -1):
    b = np.array(x_vec[((i + 1) * step) : (i + 2) * step])
    a = row_blocks[i][1] @ b
    for e in range(0,len(new_results[i])):
        x_vec[matrix_size - (num * step) + e] = new_results[i][e] - a[e]
    num = num + 1

print(x_vec)

b = np.array(results)
a = np.array(matrix)

print(np.linalg.inv(a) @ b)
            