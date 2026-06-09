import math
import re
from fractions import Fraction
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linprog


def var_key(v):
    return (0 if v.startswith('x') else 1, int(v[1:]))

def format_eq(name, const, coeffs, non_basic):
    res = ""
    if const != 0:
        res += str(const)
        
    for v in sorted(non_basic, key=var_key):
        c = coeffs.get(v, Fraction(0))
        if c == 0: continue
        
        abs_c = abs(c)
        sign = "+" if c > 0 else "-"
        c_str = "" if abs_c == 1 else str(abs_c)
        
        if res == "":
            res += f"-{c_str}{v}" if c < 0 else f"{c_str}{v}"
        else:
            res += f" {sign} {c_str}{v}"
            
    if res == "": res = "0"
    return f"{name:2} = {res}"

def parse_input():
    print("\n=== NHẬP DỮ LIỆU BÀI TOÁN TỔNG QUÁT ===")
    prob_type = input("Dạng bài toán (max / min): ").strip().lower()
    while prob_type not in ['max', 'min']:
        prob_type = input("Vui lòng nhập lại (max hoặc min): ").strip().lower()
        
    print("\nNhập hệ số hàm mục tiêu Z (ví dụ: 4x1 + 5x2 nhập là: 4 5):")
    z_inputs = list(map(Fraction, input("Hệ số: ").split()))
    num_vars = len(z_inputs)
    
    z_coeffs = {}
    for i, val in enumerate(z_inputs, 1):
        if prob_type == 'max':
            z_coeffs[f'x{i}'] = -val  
        else:
            z_coeffs[f'x{i}'] = val   
            
    num_constraints = int(input("\nNhập số lượng ràng buộc (không tính x_i >= 0): "))
    
    eqs = {'z': (Fraction(0), z_coeffs)}
    non_basic = [f'x{i}' for i in range(1, num_vars + 1)]
    basic = []
    
    geo_constraints = []
    
    print("\nNhập các ràng buộc (Ví dụ: 2x1 + 2x2 <= 9 hoặc 1x1 + 3x2 >= 5)")
    for i in range(1, num_constraints + 1):
        line = input(f"Ràng buộc {i}: ").strip()
        parts = re.split(r'\s*(<=|>=|=)\s*', line)
        
        coeffs_list = list(map(Fraction, parts[0].split()))
        sign = parts[1]
        rhs = Fraction(parts[2])
        
        geo_constraints.append({'coeffs': coeffs_list, 'sign': sign, 'rhs': rhs})
        
        if sign == '>=':
            coeffs_list = [-c for c in coeffs_list]
            rhs = -rhs
            
        w_name = f'w{i}'
        basic.append(w_name)
        
        w_coeffs = {}
        for j, val in enumerate(coeffs_list, 1):
            if val != 0:
                w_coeffs[f'x{j}'] = -val  
                
        eqs[w_name] = (rhs, w_coeffs)
        
    return eqs, basic, non_basic, prob_type, num_vars, num_constraints, geo_constraints

def display_standard_form(geo_constraints, eqs, num_vars, prob_type):
    """Hàm in ra Dạng Chuẩn của bài toán trước khi giải"""
    print(f"\n{'-'*60}")
    print("BƯỚC CHUẨN HÓA: ĐƯA BÀI TOÁN VỀ DẠNG CHUẨN (STANDARD FORM)")
    print(f"{'-'*60}")
    
    # Hàm mục tiêu
    z_str = ""
    for i in range(1, num_vars + 1):
        val = eqs['z'][1].get(f'x{i}', Fraction(0))
        val = -val if prob_type == 'max' else val
        if val == 0: continue
        sign = "+" if val > 0 and z_str != "" else ("-" if val < 0 else "")
        abs_val = abs(val)
        val_str = "" if abs_val == 1 else str(abs_val)
        z_str += f" {sign} {val_str}x{i}"
    
    z_str = z_str.strip() if z_str else "0"
    print(f"Hàm mục tiêu: Z = {z_str} -> {prob_type.upper()}")
    print("Các ràng buộc (đã thêm biến bù/thặng dư):")
    
    # Ràng buộc
    var_index = num_vars + 1
    all_vars = [f"x{i}" for i in range(1, num_vars + 1)]
    
    for idx, constr in enumerate(geo_constraints, 1):
        eq_str = ""
        coeffs = constr['coeffs']
        sign = constr['sign']
        rhs = constr['rhs']
        
        # Đảm bảo vế phải >= 0
        if rhs < 0:
            coeffs = [-c for c in coeffs]
            rhs = -rhs
            if sign == '<=': sign = '>='
            elif sign == '>=': sign = '<='
            
        for i, val in enumerate(coeffs, 1):
            if val == 0: continue
            op = "+" if val > 0 and eq_str != "" else ("-" if val < 0 else "")
            abs_val = abs(val)
            val_str = "" if abs_val == 1 else str(abs_val)
            eq_str += f" {op} {val_str}x{i}"
            
        eq_str = eq_str.strip()
        
        if sign == '<=':
            new_var = f"x{var_index}"
            eq_str += f" + {new_var}"
            all_vars.append(new_var)
            var_index += 1
        elif sign == '>=':
            new_var = f"x{var_index}"
            eq_str += f" - {new_var}"
            all_vars.append(new_var)
            var_index += 1
            
        print(f"  ({idx}) {eq_str} = {rhs}")
        
    print(f"Điều kiện không âm: {', '.join(all_vars)} >= 0")
    print(f"{'-'*60}")

def convert_data_for_other_methods(geo_constraints, num_vars, prob_type, eqs):
    z_tuple = eqs['z'][1]
    c_list = []
    for i in range(1, num_vars + 1):
        val = z_tuple.get(f'x{i}', Fraction(0))
        if prob_type == 'max':
            c_list.append(float(-val))
        else:
            c_list.append(float(val))
            
    A_ub = []
    b_ub = []
    A_eq = []
    b_eq = []
    
    for constr in geo_constraints:
        row = [float(c) for c in constr['coeffs']]
        rhs_val = float(constr['rhs'])
        
        if constr['sign'] == '<=':
            A_ub.append(row)
            b_ub.append(rhs_val)
        elif constr['sign'] == '>=':
            A_ub.append([-c for c in row])
            b_ub.append(-rhs_val)
        elif constr['sign'] == '=':
            A_eq.append(row)
            b_eq.append(rhs_val)
            
    return np.array(c_list), np.array(A_ub) if A_ub else None, np.array(b_ub) if b_ub else None, np.array(A_eq) if A_eq else None, np.array(b_eq) if b_eq else None


def solve_dictionary(eqs_orig, basic_orig, non_basic_orig, prob_type, method="Bland"):
    import copy
    eqs = copy.deepcopy(eqs_orig)
    basic = copy.deepcopy(basic_orig)
    non_basic = copy.deepcopy(non_basic_orig)
    
    print(f"\n{'='*60}")
    print(f"PHƯƠNG PHÁP: TỪ ĐIỂN ĐƠN HÌNH - {method.upper()} ({prob_type.upper()})")
    print(f"Lưu ý: Từ điển giả định điểm xuất phát cơ sở x=0 khả thi.")
    print(f"{'='*60}")
    
    iteration = 0
    while True:
        print(f"\n--- TỪ ĐIỂN {iteration} ---")
        print(format_eq('z', eqs['z'][0], eqs['z'][1], non_basic))
        for b in sorted(basic, key=var_key):
            print(format_eq(b, eqs[b][0], eqs[b][1], non_basic))
            
        z_coeffs = eqs['z'][1]
        entering = None
        
        if method == "Bland":
            for v in sorted(non_basic, key=var_key):
                if z_coeffs.get(v, 0) < 0:
                    entering = v
                    break
        else: 
            min_val = 0
            for v in non_basic:
                if z_coeffs.get(v, 0) < min_val:
                    min_val = z_coeffs[v]
                    entering = v
                    
        if entering is None:
            print("\n>> Hàm mục tiêu không còn hệ số âm. Đã đạt tối ưu!")
            break
            
        leaving = None
        min_ratio = float('inf')
        
        for b in sorted(basic, key=var_key):
            coeff = eqs[b][1].get(entering, 0)
            if coeff < 0: 
                ratio = eqs[b][0] / abs(coeff)
                if ratio < min_ratio:
                    min_ratio = ratio
                    leaving = b
                elif ratio == min_ratio and leaving is not None:
                    if var_key(b) < var_key(leaving):
                        leaving = b
                    
        if leaving is None:
            print("\n>> Bài toán không giới nội ! Z -> +vô cùng")
            return
            
        print(f"\n=> Biến XOAY: Biến VÀO = {entering}, Biến RA = {leaving}")
        
        p = eqs[leaving][1][entering]
        new_l_const = eqs[leaving][0] / -p
        
        new_l_coeffs = {leaving: Fraction(1) / p}
        for v in non_basic:
            if v != entering:
                c = eqs[leaving][1].get(v, 0)
                if c != 0: new_l_coeffs[v] = c / -p
                
        new_eqs = {}
        for k in eqs:
            if k == leaving: continue
            
            c_k, coeffs_k = eqs[k]
            a_ke = coeffs_k.get(entering, 0)
            
            new_c_k = c_k + a_ke * new_l_const
            new_coeffs_k = {}
            if a_ke != 0:
                new_coeffs_k[leaving] = a_ke / p
                
            for v in non_basic:
                if v != entering:
                    old_c = coeffs_k.get(v, 0)
                    new_c = old_c + a_ke * new_l_coeffs.get(v, 0)
                    if new_c != 0: new_coeffs_k[v] = new_c
                    
            new_eqs[k] = (new_c_k, new_coeffs_k)
            
        new_eqs[entering] = (new_l_const, new_l_coeffs)
        eqs = new_eqs
        
        basic.remove(leaving)
        basic.append(entering)
        non_basic.remove(entering)
        non_basic.append(leaving)
        
        iteration += 1

    print(f"\n{'*'*15} KẾT QUẢ TỐI ƯU ({method.upper()}) {'*'*15}")
    z_value = eqs['z'][0]
    if prob_type == 'max':
        print(f"Giá trị cực đại Z (max) = {-z_value}")
    else:
        print(f"Giá trị cực tiểu Z (min) = {z_value}")
    
    print("Nghiệm tối ưu:")
    all_vars = sorted(list(set(basic + non_basic)), key=var_key)
    for v in all_vars:
        if v.startswith('x'):
            val = eqs[v][0] if v in basic else 0
            print(f"  {v} = {val}")


def solve_scipy_2phase(c, A_ub, b_ub, A_eq, b_eq, prob_type):
    print(f"\n{'='*60}")
    print(f"PHƯƠNG PHÁP: 2 PHA (SỬ DỤNG SCIPY OPTIMIZE - HIGHS)")
    print(f"{'='*60}")
    
    c_input = c.copy() if prob_type == 'min' else -c.copy()
    x_bounds = (0, None)
    bounds = [x_bounds] * len(c)
    
    res = linprog(c_input, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
    
    if res.success:
        print(f"Trạng thái thành công: {res.message}")
        for i, val in enumerate(res.x, 1):
            print(f"Nghiệm tối ưu x{i} = {val:.4f}")
        opt_val = res.fun if prob_type == 'min' else -res.fun
        print(f"Giá trị tối ưu Z = {opt_val:.4f}")
    else:
        print(f"Không tìm thấy nghiệm tối ưu. Lý do: {res.message}")

def check_feasibility(pt, constraints, tol=1e-7):
    for constr in constraints:
        val = np.dot(constr['a'], pt)
        b = constr['b']
        if constr['sign'] == '<=' and val > b + tol: return False
        if constr['sign'] == '>=' and val < b - tol: return False
        if constr['sign'] == '=' and abs(val - b) > tol: return False
    return True

def sort_clockwise(pts):
    pts = np.unique(np.round(pts, 5), axis=0)
    if len(pts) <= 2: return pts
    centroid = np.mean(pts, axis=0)
    angles = np.arctan2(pts[:,1] - centroid[1], pts[:,0] - centroid[0])
    return pts[np.argsort(angles)]

def draw_objective_line(ax, c, z, x_lim, color, linestyle, lw, label=None):
    if abs(c[1]) > 1e-7:
        x_vals = np.array(x_lim)
        y_vals = (z - c[0] * x_vals) / c[1]
        ax.plot(x_vals, y_vals, color=color, linestyle=linestyle, linewidth=lw, label=label)
    else:
        ax.axvline(x=z/c[0], color=color, linestyle=linestyle, linewidth=lw, label=label)

def get_anchor_point(a1, a2, b_val, x_lim, y_lim):
    pts = []
    if abs(a2) > 1e-7:
        for x in x_lim:
            y = (b_val - a1 * x) / a2
            if y_lim[0] - 1e-5 <= y <= y_lim[1] + 1e-5: pts.append((x, y))
    if abs(a1) > 1e-7:
        for y in y_lim:
            x = (b_val - a2 * y) / a1
            if x_lim[0] - 1e-5 <= x <= x_lim[1] + 1e-5: pts.append((x, y))
            
    unique_pts = []
    for p in pts:
        if not any(math.hypot(p[0]-up[0], p[1]-up[1]) < 1e-5 for up in unique_pts):
            unique_pts.append(p)
            
    if len(unique_pts) >= 2: return ((unique_pts[0][0] + unique_pts[1][0])/2, (unique_pts[0][1] + unique_pts[1][1])/2)
    elif len(unique_pts) == 1: return unique_pts[0]
    else: return (0, b_val/a2) if abs(a2) > 1e-7 else (b_val/a1, 0)

def solve_geometry(opt_type, c, geo_constraints):
    print(f"\n{'='*60}")
    print(f"PHƯƠNG PHÁP: HÌNH HỌC TRỰỢT HÀM MỤC TIÊU (2 ẨN)")
    print(f"{'='*60}")
    
    opt_type = opt_type.upper()
    constraints = []
    for gc in geo_constraints:
        constraints.append({'a': np.array([float(gc['coeffs'][0]), float(gc['coeffs'][1])]), 'sign': gc['sign'], 'b': float(gc['rhs'])})
    
    constraints.append({'a': np.array([1.0, 0.0]), 'sign': '>=', 'b': 0.0})
    constraints.append({'a': np.array([0.0, 1.0]), 'sign': '>=', 'b': 0.0})
    
    M = 10000 
    calc_constraints = constraints.copy()
    calc_constraints.extend([
        {'a': np.array([1, 0]), 'sign': '<=', 'b': M},
        {'a': np.array([-1, 0]), 'sign': '<=', 'b': M},
        {'a': np.array([0, 1]), 'sign': '<=', 'b': M},
        {'a': np.array([0, -1]), 'sign': '<=', 'b': M}
    ])
    
    intersections = []
    raw_intersections = []
    
    for i, j in combinations(range(len(constraints)), 2):
        try:
            pt = np.linalg.solve([constraints[i]['a'], constraints[j]['a']], [constraints[i]['b'], constraints[j]['b']])
            raw_intersections.append(pt)
        except np.linalg.LinAlgError: pass
            
    for i, j in combinations(range(len(calc_constraints)), 2):
        try:
            pt = np.linalg.solve([calc_constraints[i]['a'], calc_constraints[j]['a']], [calc_constraints[i]['b'], calc_constraints[j]['b']])
            intersections.append(pt)
        except np.linalg.LinAlgError: pass
            
    feasible_pts = [pt for pt in intersections if check_feasibility(pt, calc_constraints)]
    
    status = ''
    p1 = p2 = None 
    draw_pts = [] 
    
    if len(feasible_pts) == 0:
        status = 'INFEASIBLE'
        print("KẾT LUẬN: BÀI TOÁN VÔ NGHIỆM (Infeasible)")
        draw_pts = raw_intersections 
    else:
        feasible_pts = sort_clockwise(feasible_pts)
        z_values = np.dot(feasible_pts, c)
        best_z = np.max(z_values) if opt_type == 'MAX' else np.min(z_values)
        optimal_pts = feasible_pts[np.abs(z_values - best_z) < 1e-5]
        is_unbounded = any(abs(abs(pt[i]) - M) < 1e-5 for pt in optimal_pts for i in range(2))
        
        if is_unbounded:
            status = 'UNBOUNDED'
            inf_val = "+ Vô cùng (+∞)" if opt_type == 'MAX' else "- Vô cùng (-∞)"
            print(f"KẾT LUẬN: BÀI TOÁN KHÔNG GIỚI NỘI (Unbounded) -> Z -> {inf_val}")
            draw_pts = [pt for pt in feasible_pts if not any(abs(abs(pt[i]) - M) < 1e-5 for i in range(2))]
        elif len(optimal_pts) > 1:
            max_dist = 0
            p1, p2 = optimal_pts[0], optimal_pts[0]
            for i, j in combinations(range(len(optimal_pts)), 2):
                dist = math.hypot(optimal_pts[i][0]-optimal_pts[j][0], optimal_pts[i][1]-optimal_pts[j][1])
                if dist > max_dist:
                    max_dist = dist
                    p1, p2 = optimal_pts[i], optimal_pts[j]
                    
            if max_dist > 1e-5:
                status = 'MULTIPLE_OPTIMAL'
                print(f"KẾT LUẬN: BÀI TOÁN VÔ SỐ NGHIỆM. Giá trị {opt_type} z* = {best_z:.2f}")
                print(f">> Đoạn nối từ A({p1[0]:.2f}, {p1[1]:.2f}) đến B({p2[0]:.2f}, {p2[1]:.2f})")
                draw_pts = feasible_pts
                best_v = p1 
            else:
                status = 'OPTIMAL'
                best_v = optimal_pts[0]
        else: 
            status = 'OPTIMAL'
            best_v = optimal_pts[0]
                
        if status == 'OPTIMAL':
            print("KẾT LUẬN: BÀI TOÁN CÓ NGHIỆM DUY NHẤT")
            print(f">> Nghiệm: x1 = {best_v[0]:.2f}, x2 = {best_v[1]:.2f}")
            print(f">> Giá trị cực trị {opt_type} z = {best_z:.2f}")
            draw_pts = feasible_pts

    fig, ax = plt.subplots(figsize=(9, 7))
    draw_pts = np.array(draw_pts)
    
    if len(draw_pts) > 0:
        x_min, x_max = np.min(draw_pts[:,0]), np.max(draw_pts[:,0])
        y_min, y_max = np.min(draw_pts[:,1]), np.max(draw_pts[:,1])
        margin_x, margin_y = max(2, (x_max - x_min)*0.4), max(2, (y_max - y_min)*0.4)
    else:
        x_min, x_max, y_min, y_max = 0, 10, 0, 10
        margin_x = margin_y = 2

    x_lim_draw, y_lim_draw = [x_min - margin_x, x_max + margin_x], [y_min - margin_y, y_max + margin_y]
    arrow_len = (x_max - x_min + 2*margin_x) * 0.08 

    for idx, constr in enumerate(constraints):
        a1, a2, b_val = constr['a'][0], constr['a'][1], constr['b']
        if abs(a2) > 1e-7:
            x_vals = np.array(x_lim_draw)
            ax.plot(x_vals, (b_val - a1 * x_vals) / a2, color='black', lw=1.5, alpha=0.6)
        else:
            ax.axvline(x=b_val/a1, color='black', lw=1.5, alpha=0.6)
            
        anchor = get_anchor_point(a1, a2, b_val, x_lim_draw, y_lim_draw)
        ax.text(anchor[0], anchor[1], f'({idx+1})', color='black', fontsize=10, fontweight='bold', 
                ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', boxstyle='circle'))

    if status != 'INFEASIBLE' and len(feasible_pts) > 2:
        polygon = plt.Polygon(feasible_pts, closed=True, color='cyan', alpha=0.3, label='Mien nghiem (Feasible region)')
        ax.add_patch(polygon)
    
    zx, zy = c[0], c[1]
    if opt_type == 'MIN': zx, zy = -zx, -zy 
    norm_z = math.hypot(zx, zy)
    if norm_z > 1e-7: zx, zy = (zx / norm_z) * arrow_len * 1.5, (zy / norm_z) * arrow_len * 1.5

    if status in ['OPTIMAL', 'MULTIPLE_OPTIMAL']:
        worst_idx = np.argmin(z_values) if opt_type == 'MAX' else np.argmax(z_values)
        worst_z = z_values[worst_idx]
        for z_val in np.linspace(worst_z, best_z, 4)[:-1]:
            draw_objective_line(ax, c, z_val, x_lim_draw, color='gray', linestyle=':', lw=1.5)
        draw_objective_line(ax, c, best_z, x_lim_draw, color='red', linestyle='--', lw=2.5, label=f'Duong muc z*={best_z:.2f}')
        
        if status == 'MULTIPLE_OPTIMAL':
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='magenta', lw=5, zorder=6, label='Tap vo so nghiem')
            vector_origin = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2)
        else:
            ax.plot(best_v[0], best_v[1], marker='*', color='red', markersize=15, zorder=5)
            vector_origin = (best_v[0], best_v[1])
        
        ax.annotate('', xy=(vector_origin[0] + zx, vector_origin[1] + zy), xytext=(vector_origin[0], vector_origin[1]),
                    arrowprops=dict(arrowstyle="simple,head_width=0.6,head_length=0.8", color='red', lw=1.5))

    ax.set_xlim(x_lim_draw)
    ax.set_ylim(y_lim_draw)
    ax.axhline(0, color='black', lw=1)
    ax.axvline(0, color='black', lw=1)
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_title(f'Khao sat hinh hoc: {opt_type} Z', fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    if status in ['OPTIMAL', 'MULTIPLE_OPTIMAL']: ax.legend()
    plt.tight_layout()
    print(">> Đang hiển thị biểu đồ trực quan hình học...")
    plt.show()


def main():
    try:
        eqs, basic, non_basic, prob_type, num_vars, num_constraints, geo_constraints = parse_input()
        c_np, A_ub, b_ub, A_eq, b_eq = convert_data_for_other_methods(geo_constraints, num_vars, prob_type, eqs)
        
        while True:
            print("\n" + "="*45)
            print("DANH SÁCH PHƯƠNG PHÁP GIẢI QUY HOẠCH TUYẾN TÍNH:")
            print("1. Phương pháp đơn hình Từ điển (Thuật toán Dantzig)")
            print("2. Phương pháp đơn hình Từ điển chống vòng lặp (Quy tắc Bland)")
            print("3. Phương pháp 2 Pha / Gọi thư viện Scipy HiGHS")
            print("4. Phương pháp hình học (Chỉ khả thi cho bài toán 2 ẩn)")
            print("5. Giải bằng TẤT CẢ 4 CÁCH trên")
            print("6. Thoát chương trình")
            print("="*45)
            
            choice = input("Chọn cách giải của bạn (1->6): ").strip()
            
            if choice in ['1', '2', '3', '5']:
                # Hiển thị dạng chuẩn trước khi giải các thuật toán đại số
                display_standard_form(geo_constraints, eqs, num_vars, prob_type)
            
            if choice == '1':
                solve_dictionary(eqs, basic, non_basic, prob_type, method="Don hinh")
            elif choice == '2':
                solve_dictionary(eqs, basic, non_basic, prob_type, method="Bland")
            elif choice == '3':
                solve_scipy_2phase(c_np, A_ub, b_ub, A_eq, b_eq, prob_type)
            elif choice == '4':
                if num_vars != 2:
                    print(f"\n[Lỗi] Phương pháp hình học chỉ áp dụng cho bài toán có 2 ẩn số. Bài toán hiện tại có {num_vars} ẩn.")
                else:
                    solve_geometry(prob_type, c_np, geo_constraints)
            elif choice == '5':
                print("\n>>> BẮT ĐẦU CHẠY THỬ NGHIỆM ĐỒNG THỜI 4 PHƯƠNG PHÁP <<<")
                solve_dictionary(eqs, basic, non_basic, prob_type, method="Don hinh")
                solve_dictionary(eqs, basic, non_basic, prob_type, method="Bland")
                solve_scipy_2phase(c_np, A_ub, b_ub, A_eq, b_eq, prob_type)
                if num_vars == 2:
                    solve_geometry(prob_type, c_np, geo_constraints)
                else:
                    print(f"\n[Bỏ qua hình học] Số biến = {num_vars} (!= 2)")
            elif choice == '6':
                print("\nCảm ơn bạn đã sử dụng chương trình!")
                break
            else:
                print("Lựa chọn không hợp lệ, vui lòng chọn lại từ 1 đến 6.")
                
    except Exception as e:
        print(f"\n[Hệ thống gặp lỗi]: {e}")
        print("Vui lòng kiểm tra lại cấu trúc cú pháp nhập vào.")

if __name__ == "__main__":
    main()