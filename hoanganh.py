from scipy.optimize import linprog
c = [-1, -3]
A = [
    [-1, -1],
    [-1, 1],
    [1, 2]
]
B = [-3, -1, 4]
x_bounds = (0, None)
res = linprog(c, A_ub=A, b_ub=B, bounds=[x_bounds, x_bounds], method='highs')
if res.success:
    print(f"Trạng thái: {res.message}")
    print(f"Nghiệm tối ưu x1 = {res.x[0]}")
    print(f"Nghiệm tối ưu x2 = {res.x[1]}")
    print(f"Giá trị tối ưu z = {res.fun}")
else:
    print("Không tìm thấy nghiệm tối ưu.")