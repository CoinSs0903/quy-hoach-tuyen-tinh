import numpy as np
import sys
from typing import List, Dict, Tuple, Optional

def get_user_input() -> Optional[Tuple[str, int, np.ndarray, List[Dict]]]:
    """Hàm thu thập và kiểm tra dữ liệu đầu vào cho bài toán n biến."""
    print("=== GIẢI QHTT: THUẬT TOÁN ĐƠN HÌNH - QUY TẮC BLAND ===")
    
    opt_type = input("Loại bài toán (MAX hoặc MIN): ").strip().upper()
    if opt_type not in ['MAX', 'MIN']:
        print("Lỗi: Chỉ chấp nhận MAX hoặc MIN.")
        return None
        
    try:
        n = int(input("Nhập số lượng biến quyết định n: "))
        if n < 1:
            raise ValueError("Số biến phải lớn hơn 0.")
            
        c_input = input(f"Nhập hệ số hàm mục tiêu c (gồm {n} số, cách nhau khoảng trắng): ")
        c = np.array([float(val) for val in c_input.split()])
        if len(c) != n:
            raise ValueError(f"Phải nhập chính xác {n} hệ số.")
    except Exception as e:
        print(f"Lỗi nhập hàm mục tiêu: {e}")
        return None
    
    try:
        m = int(input("Nhập số lượng ràng buộc m (không tính x >= 0): "))
        constraints_info = []
        
        print(f"\nNhập ràng buộc theo dạng: a1 a2 ... an dấu b (VD: 2 1 -1 <= 10)")
        for i in range(m):
            constr = input(f"Ràng buộc {i+1}: ").strip().split()
            if len(constr) != n + 2:
                print(f"Lỗi: Ràng buộc phải gồm {n} hệ số, 1 dấu và 1 hệ số tự do b.")
                return None
                
            a_coeffs = [float(x) for x in constr[:n]]
            sign = constr[n]
            b_val = float(constr[n+1])
            
            if sign not in ['<=', '>=', '=']:
                print("Lỗi: Dấu ràng buộc phải là '<=', '>=', hoặc '='.")
                return None
                
            constraints_info.append({'a': np.array(a_coeffs), 'sign': sign, 'b': b_val})
    except Exception as e:
        print(f"Lỗi nhập ràng buộc: {e}")
        return None
        
    return opt_type, n, c, constraints_info

def print_tableau(T: np.ndarray, basis: List[int], N: int, m: int, iteration: int):
    """Hàm in bảng Đơn hình (Tableau) ra màn hình với định dạng cột."""
    if iteration == 0:
        print("\n[BẢNG ĐƠN HÌNH BAN ĐẦU]")
    else:
        print(f"\n[BẢNG ĐƠN HÌNH SAU VÒNG LẶP {iteration}]")

    # In tiêu đề cột
    header = f"{'Cơ sở':>7} |"
    for i in range(N):
        header += f" x{i+1:>6} |"
    header += f" {'RHS':>8}"
    
    print("-" * len(header))
    print(header)
    print("-" * len(header))

    # In dòng hàm mục tiêu (Row 0)
    row0 = f"{'z (D0)':>7} |"
    for val in T[0, :N]:
        row0 += f" {val:>7.2f} |"
    row0 += f" {T[0, N]:>8.2f}"
    print(row0)
    print("-" * len(header))

    # In các dòng ràng buộc (Row 1 đến m)
    for i in range(1, m + 1):
        base_var = f"x{basis[i-1]+1}" if basis[i-1] != -1 else "Art"
        row_str = f"{base_var} (D{i}) |"
        for val in T[i, :N]:
            row_str += f" {val:>7.2f} |"
        row_str += f" {T[i, N]:>8.2f}"
        print(row_str)
        
    print("-" * len(header))

def solve_simplex_bland(opt_type: str, n: int, c: np.ndarray, constraints: List[Dict]):
    """Giải bài toán QHTT bằng thuật toán Đơn hình, in ra từng bước xoay."""
    m = len(constraints)
    var_count = n
    artificial_cols = []
    
    for constr in constraints:
        if constr['sign'] == '<=': var_count += 1
        elif constr['sign'] == '>=': var_count += 2
        elif constr['sign'] == '=': var_count += 1
            
    N = var_count 
    T = np.zeros((m + 1, N + 1))
    basis = [-1] * m
    
    c_max = c if opt_type == 'MAX' else -c
    T[0, 0:n] = -c_max 
    col_idx = n
    M = 10000 
    
    for i, constr in enumerate(constraints):
        a = constr['a'].copy()
        b = constr['b']
        sign = constr['sign']
        
        if b < 0:
            a, b = -a, -b
            if sign == '<=': sign = '>='
            elif sign == '>=': sign = '<='
            
        T[i+1, 0:n] = a
        T[i+1, N] = b
        
        if sign == '<=':
            T[i+1, col_idx] = 1
            basis[i] = col_idx
            col_idx += 1
        elif sign == '>=':
            T[i+1, col_idx] = -1      
            T[i+1, col_idx+1] = 1    
            T[0, col_idx+1] = M      
            basis[i] = col_idx + 1
            artificial_cols.append(col_idx + 1)
            col_idx += 2
        elif sign == '=':
            T[i+1, col_idx] = 1      
            T[0, col_idx] = M        
            basis[i] = col_idx
            artificial_cols.append(col_idx)
            col_idx += 1
            
    # Hiệu chỉnh Row 0 (Khử M)
    for i in range(m):
        if basis[i] in artificial_cols:
            T[0, :] -= M * T[i+1, :]
            
    # IN BẢNG BAN ĐẦU
    print_tableau(T, basis, N, m, 0)
            
    # VÒNG LẶP ĐƠN HÌNH
    iteration = 0
    status = ""
    while True:
        iteration += 1
        
        # BƯỚC 1: Chọn biến vào
        enter_col = -1
        for j in range(N):
            if T[0, j] < -1e-7:
                enter_col = j
                break 
                
        if enter_col == -1:
            status = "OPTIMAL"
            break
            
        # BƯỚC 2: Chọn biến ra
        min_ratio = float('inf')
        leave_row = -1
        leave_var_index = float('inf') 
        
        for i in range(1, m + 1):
            if T[i, enter_col] > 1e-7:
                ratio = T[i, N] / T[i, enter_col]
                if ratio < min_ratio - 1e-7:
                    min_ratio = ratio
                    leave_row = i
                    leave_var_index = basis[i-1]
                elif abs(ratio - min_ratio) <= 1e-7:
                    if basis[i-1] < leave_var_index:
                        leave_row = i
                        leave_var_index = basis[i-1]
                        
        if leave_row == -1:
            status = "UNBOUNDED"
            break
            
        # --- IN THÔNG TIN BƯỚC XOAY VÀ PHÉP TOÁN ---
        print(f"\n=> ĐANG XÉT VÒNG LẶP {iteration}:")
        print(f"   + Cột xoay (Biến vào): x{enter_col + 1}")
        print(f"   + Dòng xoay (Biến ra): x{basis[leave_row - 1] + 1} (Dòng {leave_row})")
        
        pivot_val = T[leave_row, enter_col]
        print(f"   + Phần tử trục (Pivot): {pivot_val:.2f}")
        print("\n   [CÁC PHÉP BIẾN ĐỔI DÒNG GAUSS-JORDAN]")
        
        if abs(pivot_val - 1.0) > 1e-7:
            print(f"   -> D{leave_row} (mới) = D{leave_row} (cũ) / {pivot_val:.2f}")
        else:
            print(f"   -> D{leave_row} (mới) giữ nguyên do Pivot = 1.00")
            
        # BƯỚC 3: Xoay bảng
        T[leave_row, :] /= pivot_val
        for i in range(m + 1):
            if i != leave_row:
                factor = T[i, enter_col]
                if abs(factor) > 1e-7:
                    sign_str = "-" if factor > 0 else "+"
                    print(f"   -> D{i} (mới) = D{i} (cũ) {sign_str} {abs(factor):.2f} * D{leave_row} (mới)")
                    T[i, :] -= factor * T[leave_row, :]
                
        # Cập nhật cơ sở và in bảng mới
        basis[leave_row - 1] = enter_col
        print_tableau(T, basis, N, m, iteration)

    # KẾT LUẬN
    print("\n" + "="*50)
    print(f"TỔNG SỐ VÒNG LẶP: {iteration - 1}")
    
    if status == "OPTIMAL":
        for i in range(m):
            if basis[i] in artificial_cols and T[i+1, N] > 1e-7:
                status = "INFEASIBLE"
                break

    if status == "INFEASIBLE":
        print("KẾT LUẬN: BÀI TOÁN VÔ NGHIỆM (Infeasible)")
        print(">> Các ràng buộc mâu thuẫn, biến giả vẫn nằm trong cơ sở có giá trị > 0.")
    elif status == "UNBOUNDED":
        print("KẾT LUẬN: BÀI TOÁN KHÔNG GIỚI NỘI (Unbounded)")
        inf_val = "+∞" if opt_type == 'MAX' else "-∞"
        print(f">> Hàm mục tiêu có thể tiến tới {inf_val}.")
    elif status == "OPTIMAL":
        print("KẾT LUẬN: ĐÃ TÌM THẤY PHƯƠNG ÁN TỐI ƯU")
        x_opt = [0.0] * n
        for i in range(m):
            if basis[i] < n: x_opt[basis[i]] = T[i+1, N]
                
        for i in range(n):
            print(f">> x_{i+1} = {x_opt[i]:.4f}")
            
        best_z = T[0, N] if opt_type == 'MAX' else -T[0, N]
        print(f">> Giá trị {opt_type} z* = {best_z:.4f}")
    print("="*50)

def main():
    inputs = get_user_input()
    if inputs:
        solve_simplex_bland(*inputs)

if __name__ == "__main__":
    main()