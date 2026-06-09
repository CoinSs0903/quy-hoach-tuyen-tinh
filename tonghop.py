import math
import re
from fractions import Fraction
from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linprog


def var_key(v):
    if v in ['z', 'z_aux']:
        return (-1, 0, '')
    type_prefix = 0 if v.startswith('x') else 1
    match = re.match(r'^([xw])(\d+)(.*)$', v)
    if match:
        num = int(match.group(2))
        suffix = match.group(3)
        return (type_prefix, num, suffix)
    return (type_prefix, 999, v)

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
    
    print("\nNhập điều kiện dấu cho từng biến số x_i (nhập '>=0' hoặc '<=0' hoặc 'td' cho tự do):")
    var_signs = []
    for i in range(1, num_vars + 1):
        sign = input(f"Điều kiện của x{i} (Mặc định >=0): ").strip().lower()
        if sign == '<=0':
            var_signs.append('<=0')
        elif sign in ['td', 'tu do', 'tự do', 'free']:
            var_signs.append('free')
        else:
            var_signs.append('>=0')
            
    # Build standard variables mapping
    var_mapping = {}
    for i in range(1, num_vars + 1):
        sign = var_signs[i - 1]
        if sign == '>=0':
            var_mapping[f'x{i}'] = [(f'x{i}', Fraction(1))]
        elif sign == '<=0':
            var_mapping[f'x{i}'] = [(f"x{i}'", Fraction(-1))]
        elif sign == 'free':
            var_mapping[f'x{i}'] = [(f'x{i}+', Fraction(1)), (f'x{i}-', Fraction(-1))]

    z_coeffs = {}
    non_basic = []
    for i, val in enumerate(z_inputs, 1):
        if prob_type == 'max':
            base_val = -val
        else:
            base_val = val
            
        mapping = var_mapping[f'x{i}']
        for std_var, factor in mapping:
            z_coeffs[std_var] = z_coeffs.get(std_var, Fraction(0)) + base_val * factor
            if std_var not in non_basic:
                non_basic.append(std_var)
                
    non_basic.sort(key=var_key)
            
    num_constraints = int(input("\nNhập số lượng ràng buộc (không tính x_i >= 0): "))
    
    geo_constraints = []
    
    print("\nNhập các ràng buộc (Ví dụ: 2x1 + 2x2 <= 9 hoặc 1x1 + 3x2 >= 5)")
    flat_constraints_data = []
    for i in range(1, num_constraints + 1):
        line = input(f"Ràng buộc {i}: ").strip()
        parts = re.split(r'\s*(<=|>=|=)\s*', line)
        
        coeffs_list = list(map(Fraction, parts[0].split()))
        sign = parts[1]
        rhs = Fraction(parts[2])
        
        geo_constraints.append({'coeffs': coeffs_list, 'sign': sign, 'rhs': rhs})
        
        if sign == '=':
            flat_constraints_data.append({'coeffs': coeffs_list, 'sign': '<=', 'rhs': rhs})
            flat_constraints_data.append({'coeffs': coeffs_list, 'sign': '>=', 'rhs': rhs})
        else:
            flat_constraints_data.append({'coeffs': coeffs_list, 'sign': sign, 'rhs': rhs})
            
    eqs = {'z': (Fraction(0), z_coeffs)}
    basic = []
    
    for i, constr in enumerate(flat_constraints_data, 1):
        coeffs_list = constr['coeffs']
        sign = constr['sign']
        rhs = constr['rhs']
        
        # Build standard coefficients in terms of std_vars
        std_coeffs = {}
        for j, val in enumerate(coeffs_list, 1):
            if val == 0: continue
            mapping = var_mapping[f'x{j}']
            for std_var, factor in mapping:
                std_coeffs[std_var] = std_coeffs.get(std_var, Fraction(0)) + val * factor
                
        factor = Fraction(1)
        if sign == '>=':
            factor = Fraction(-1)
            rhs = -rhs
            
        w_name = f'w{i}'
        basic.append(w_name)
        
        w_coeffs = {}
        for std_var, val in std_coeffs.items():
            final_val = val * factor
            if final_val != 0:
                w_coeffs[std_var] = -final_val
                
        eqs[w_name] = (rhs, w_coeffs)
        
    return eqs, basic, non_basic, prob_type, num_vars, num_constraints, geo_constraints, var_signs, z_inputs

def display_standard_form(geo_constraints, eqs, num_vars, prob_type):
    """Hàm in ra Dạng Chuẩn của bài toán trước khi giải"""
    print(f"\n{'-'*60}")
    print("BƯỚC CHUẨN HÓA: ĐƯA BÀI TOÁN VỀ DẠNG CHUẨN (STANDARD FORM)")
    print(f"{'-'*60}")
    
    # Hàm mục tiêu
    z_str = ""
    std_vars = sorted(list(eqs['z'][1].keys()), key=var_key)
    for v in std_vars:
        val = eqs['z'][1].get(v, Fraction(0))
        val = -val if prob_type == 'max' else val
        if val == 0: continue
        sign = "+" if val > 0 and z_str != "" else ("-" if val < 0 else "")
        abs_val = abs(val)
        val_str = "" if abs_val == 1 else str(abs_val)
        z_str += f" {sign} {val_str}{v}"
        
    z_str = z_str.strip() if z_str else "0"
    print(f"Hàm mục tiêu: Z = {z_str} -> {prob_type.upper()}")
    print("Các ràng buộc dạng chuẩn (biến chuẩn hóa >= 0):")
    
    # Ràng buộc
    all_vars = sorted(list(eqs['z'][1].keys()), key=var_key)
    idx = 1
    for b in sorted(eqs.keys(), key=var_key):
        if b in ['z', 'z_aux']: continue
        const, coeffs = eqs[b]
        row_terms = []
        for v in sorted(coeffs.keys(), key=var_key):
            val = -coeffs[v]
            if val == 0: continue
            sign = "+" if val > 0 and len(row_terms) > 0 else ("-" if val < 0 else "")
            abs_val = abs(val)
            val_str = "" if abs_val == 1 else str(abs_val)
            row_terms.append(f"{sign} {val_str}{v}".strip())
        
        row_terms.append(f"+ {b}")
        if b not in all_vars:
            all_vars.append(b)
            
        eq_str = " ".join(row_terms).replace(" + -", " - ").replace(" + +", " + ")
        if eq_str.startswith("+ "):
            eq_str = eq_str[2:]
        print(f"  ({idx}) {eq_str} = {const}")
        idx += 1
        
    all_vars.sort(key=var_key)
    print(f"Điều kiện không âm: {', '.join(all_vars)} >= 0")
    print(f"{'-'*60}")

def convert_data_for_other_methods(geo_constraints, num_vars, prob_type, z_inputs):
    c_list = []
    for val in z_inputs:
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


def reconstruct_original_solution_py(eqs, basic, non_basic, var_signs, num_vars):
    original_solution = {}
    for i in range(1, num_vars + 1):
        sign = var_signs[i - 1] if var_signs else '>=0'
        
        if sign == '>=0':
            v_name = f'x{i}'
            val = eqs[v_name][0] if v_name in basic else Fraction(0)
        elif sign == '<=0':
            v_name = f"x{i}'"
            val = -(eqs[v_name][0] if v_name in basic else Fraction(0))
        elif sign == 'free':
            v_plus = f'x{i}+'
            v_minus = f'x{i}-'
            val_plus = eqs[v_plus][0] if v_plus in basic else Fraction(0)
            val_minus = eqs[v_minus][0] if v_minus in basic else Fraction(0)
            val = val_plus - val_minus
            
        original_solution[f'x{i}'] = val
    return original_solution

def check_multiple_optimal_py(eqs, basic, non_basic, var_signs, num_vars):
    if not var_signs or not num_vars:
        return None
    z_coeffs = eqs['z'][1]
    candidates = []
    for v in non_basic:
        if z_coeffs.get(v, Fraction(0)) == 0:
            candidates.append(v)
            
    if not candidates:
        return None
        
    for entering in candidates:
        leaving = None
        min_ratio = float('inf')
        for b in basic:
            coeff = eqs[b][1].get(entering, Fraction(0))
            if coeff < 0:
                ratio = eqs[b][0] / abs(coeff)
                if ratio < min_ratio:
                    min_ratio = ratio
                    leaving = b
                elif ratio == min_ratio and leaving is not None:
                    if var_key(b) < var_key(leaving):
                        leaving = b
                        
        if leaving is not None:
            p = eqs[leaving][1][entering]
            entering_val = eqs[leaving][0] / -p
            
            new_vals = {}
            new_vals[entering] = entering_val
            new_vals[leaving] = Fraction(0)
            
            for b in basic:
                if b != leaving:
                    coeff_ent = eqs[b][1].get(entering, Fraction(0))
                    new_vals[b] = eqs[b][0] + coeff_ent * entering_val
            for nb in non_basic:
                if nb != entering:
                    new_vals[nb] = Fraction(0)
                    
            sol2 = {}
            for i in range(1, num_vars + 1):
                sign = var_signs[i - 1]
                if sign == '>=0':
                    val = new_vals.get(f'x{i}', Fraction(0))
                elif sign == '<=0':
                    val = -new_vals.get(f"x{i}'", Fraction(0))
                elif sign == 'free':
                    val = new_vals.get(f'x{i}+', Fraction(0)) - new_vals.get(f'x{i}-', Fraction(0))
                sol2[f'x{i}'] = val
            return sol2
    return None

def is_problem_feasible_py(eqs, basic, non_basic):
    from scipy.optimize import linprog
    all_vars = sorted(list(non_basic) + list(basic), key=var_key)
    var_to_idx = {v: i for i, v in enumerate(all_vars)}
    A_eq = []
    b_eq = []
    for b in basic:
        row = [0] * len(all_vars)
        row[var_to_idx[b]] = 1
        const, coeffs = eqs[b]
        for nb, coeff in coeffs.items():
            if nb in var_to_idx:
                row[var_to_idx[nb]] = -float(coeff)
        A_eq.append(row)
        b_eq.append(float(const))
    c = [0] * len(all_vars)
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method='highs')
    return res.success, res

def solve_feasible_dictionary_py(eqs, basic, non_basic, prob_type, var_signs, num_vars):
    from scipy.optimize import linprog
    all_vars = sorted(list(non_basic) + list(basic), key=var_key)
    var_to_idx = {v: i for i, v in enumerate(all_vars)}
    A_eq = []
    b_eq = []
    for b in basic:
        row = [0] * len(all_vars)
        row[var_to_idx[b]] = 1
        const, coeffs = eqs[b]
        for nb, coeff in coeffs.items():
            if nb in var_to_idx:
                row[var_to_idx[nb]] = -float(coeff)
        A_eq.append(row)
        b_eq.append(float(const))
    c = [0] * len(all_vars)
    z_const, z_coeffs = eqs['z']
    for nb, coeff in z_coeffs.items():
        c[var_to_idx[nb]] = float(coeff) if prob_type == 'min' else -float(coeff)
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=(0, None), method='highs')
    return res


def solve_dictionary(eqs_orig, basic_orig, non_basic_orig, prob_type, method="Bland", var_signs=None, num_vars=None):
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
            print(f"\n{'='*20} KẾT LUẬN: BÀI TOÁN KHÔNG GIỚI NỘI (UNBOUNDED) {'='*20}")
            z_limit = "+∞" if prob_type == 'max' else "-∞"
            z_label = "Z (max)" if prob_type == 'max' else "Z (min)"
            print(f"Hàm mục tiêu {z_label} -> {z_limit}")
            print("Bài toán không có nghiệm tối ưu hữu hạn.")
            print(f"{'='*60}")
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

    # Check if there is any negative basic variable in standard dictionary
    has_negative_basic = any(eqs[b][0] < 0 for b in basic)
    if has_negative_basic:
        # Check if the problem is actually feasible
        is_feasible, res_feas = is_problem_feasible_py(eqs, basic, non_basic)
        if is_feasible:
            res_opt = solve_feasible_dictionary_py(eqs, basic, non_basic, prob_type, var_signs, num_vars)
            print(f"\n{'='*20} KẾT LUẬN: BÀI TOÁN CÓ NGHIỆM TỐI ƯU {'='*20}")
            opt_val = res_opt.fun + float(eqs['z'][0]) if prob_type == 'min' else -res_opt.fun + float(eqs['z'][0])
            opt_label = "cực tiểu Z (min)" if prob_type == 'min' else "cực đại Z (max)"
            print(f"Giá trị {opt_label} = {opt_val:.4f}")
            
            all_vars = sorted(list(non_basic) + list(basic), key=var_key)
            var_values = {v: res_opt.x[i] for i, v in enumerate(all_vars)}
            
            orig_sol = {}
            if var_signs and num_vars:
                for i in range(1, num_vars + 1):
                    sign = var_signs[i - 1]
                    if sign == '>=0':
                        orig_sol[f'x{i}'] = var_values.get(f'x{i}', 0.0)
                    elif sign == '<=0':
                        orig_sol[f'x{i}'] = -var_values.get(f"x{i}'", 0.0)
                    elif sign == 'free':
                        orig_sol[f'x{i}'] = var_values.get(f'x{i}+', 0.0) - var_values.get(f'x{i}-', 0.0)
            
            if orig_sol:
                var_names = ", ".join(f"x{i}" for i in range(1, num_vars + 1))
                var_vals = ", ".join(f"{orig_sol[f'x{i}']:.4f}" for i in range(1, num_vars + 1))
                print(f"Nghiệm tối ưu: x* = ({var_names}) = ({var_vals})")
            print(f"{'='*60}")
            return
        else:
            print(f"\n{'='*20} KẾT LUẬN: BÀI TOÁN VÔ NGHIỆM (INFEASIBLE) {'='*20}")
            neg_vars = ", ".join(f"{b} = {eqs[b][0]}" for b in sorted(basic, key=var_key) if eqs[b][0] < 0)
            print(f"Từ điển tối ưu nhưng không khả thi do chứa các biến cơ sở âm ({neg_vars}).")
            print("Vui lòng tham khảo phương pháp 2 Pha.")
            opt_label = "Z (max)" if prob_type == 'max' else "Z (min)"
            opt_val_disp = "-∞" if prob_type == 'max' else "+∞"
            print(f"Hàm mục tiêu {opt_label} = {opt_val_disp}")
            print(f"{'='*60}")
            return

    # Check for multiple optimal solutions
    orig_sol2 = None
    if var_signs and num_vars:
        orig_sol1 = reconstruct_original_solution_py(eqs, basic, non_basic, var_signs, num_vars)
        orig_sol2 = check_multiple_optimal_py(eqs, basic, non_basic, var_signs, num_vars)
        
        # Double check if they are actually different
        if orig_sol2:
            diff = False
            for i in range(1, num_vars + 1):
                if orig_sol1[f'x{i}'] != orig_sol2[f'x{i}']:
                    diff = True
                    break
            if not diff:
                orig_sol2 = None

    z_value = eqs['z'][0]
    opt_label = "cực đại Z (max)" if prob_type == 'max' else "cực tiểu Z (min)"
    opt_val_disp = -z_value if prob_type == 'max' else z_value

    if orig_sol2 is not None:
        print(f"\n{'='*20} KẾT LUẬN: BÀI TOÁN VÔ SỐ NGHIỆM ({method.upper()}) {'='*20}")
        print(f"Giá trị tối ưu {opt_label} = {opt_val_disp}")
        var_names = ", ".join(f"x{i}" for i in range(1, num_vars + 1))
        var_vals1 = ", ".join(str(orig_sol1[f'x{i}']) for i in range(1, num_vars + 1))
        var_vals2 = ", ".join(str(orig_sol2[f'x{i}']) for i in range(1, num_vars + 1))
        print(f"Nghiệm cực biên tối ưu thứ nhất: X1 = ({var_names}) = ({var_vals1})")
        print(f"Nghiệm cực biên tối ưu thứ hai: X2 = ({var_names}) = ({var_vals2})")
        print("Nghiệm tổng quát của bài toán (với 0 <= λ <= 1):")
        for i in range(1, num_vars + 1):
            val1 = orig_sol1[f'x{i}']
            val2 = orig_sol2[f'x{i}']
            if val1 == val2:
                print(f"  x{i} = {val1}")
            else:
                print(f"  x{i} = {val1} * λ + {val2} * (1 - λ)")
    else:
        print(f"\n{'='*20} KẾT LUẬN: BÀI TOÁN CÓ NGHIỆM TỐI ƯU ({method.upper()}) {'='*20}")
        print(f"Giá trị {opt_label} = {opt_val_disp}")
        if var_signs and num_vars:
            var_names = ", ".join(f"x{i}" for i in range(1, num_vars + 1))
            var_vals = ", ".join(str(orig_sol1[f'x{i}']) for i in range(1, num_vars + 1))
            print(f"Nghiệm tối ưu: x* = ({var_names}) = ({var_vals})")
        
    all_vars = sorted(list(set(basic + non_basic)), key=var_key)
    dict_vals = []
    for v in all_vars:
        val = eqs[v][0] if v in basic else 0
        dict_vals.append(f"{v} = {val}")
    print(f"Biến chuẩn hóa: {', '.join(dict_vals)}")
    print(f"{'='*60}")


def solve_scipy_2phase(c, A_ub, b_ub, A_eq, b_eq, prob_type, var_signs=None):
    print(f"\n{'='*60}")
    print(f"PHƯƠNG PHÁP: 2 PHA (SỬ DỤNG SCIPY OPTIMIZE - HIGHS)")
    print(f"{'='*60}")
    
    c_input = c.copy() if prob_type == 'min' else -c.copy()
    if var_signs:
        bounds = []
        for sign in var_signs:
            if sign == '>=0':
                bounds.append((0, None))
            elif sign == '<=0':
                bounds.append((None, 0))
            elif sign == 'free':
                bounds.append((None, None))
            else:
                bounds.append((0, None))
    else:
        bounds = [(0, None)] * len(c)
    
    res = linprog(c_input, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
    
    if res.success:
        print(f"\n{'='*20} KẾT LUẬN: BÀI TOÁN CÓ NGHIỆM TỐI ƯU {'='*20}")
        opt_val = res.fun if prob_type == 'min' else -res.fun
        opt_label = "cực tiểu Z (min)" if prob_type == 'min' else "cực đại Z (max)"
        print(f"Giá trị {opt_label} = {opt_val:.4f}")
        
        # Coordinate vector
        var_names = ", ".join(f"x{i}" for i in range(1, len(res.x) + 1))
        var_vals = ", ".join(f"{val:.4f}" for val in res.x)
        print(f"Nghiệm tối ưu: x* = ({var_names}) = ({var_vals})")
        print(f"{'='*60}")
    else:
        print(f"\n{'='*20} KẾT LUẬN: BÀI TOÁN VÔ NGHIỆM (INFEASIBLE) {'='*20}")
        print(f"Không tìm thấy nghiệm tối ưu. Lý do: {res.message}")
        opt_label = "Z (max)" if prob_type == 'max' else "Z (min)"
        opt_val_disp = "-∞" if prob_type == 'max' else "+∞"
        print(f"Hàm mục tiêu {opt_label} = {opt_val_disp}")
        print(f"{'='*60}")

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

def solve_geometry(opt_type, c, geo_constraints, var_signs=None):
    print(f"\n{'='*60}")
    print(f"PHƯƠNG PHÁP: HÌNH HỌC TRỰỢT HÀM MỤC TIÊU (2 ẨN)")
    print(f"{'='*60}")
    
    opt_type = opt_type.upper()
    constraints = []
    for gc in geo_constraints:
        constraints.append({'a': np.array([float(gc['coeffs'][0]), float(gc['coeffs'][1])]), 'sign': gc['sign'], 'b': float(gc['rhs'])})
    
    # Add bounds based on var_signs
    if var_signs:
        if var_signs[0] == '>=0':
            constraints.append({'a': np.array([1.0, 0.0]), 'sign': '>=', 'b': 0.0})
        elif var_signs[0] == '<=0':
            constraints.append({'a': np.array([1.0, 0.0]), 'sign': '<=', 'b': 0.0})
        # If free, do nothing
        
        if var_signs[1] == '>=0':
            constraints.append({'a': np.array([0.0, 1.0]), 'sign': '>=', 'b': 0.0})
        elif var_signs[1] == '<=0':
            constraints.append({'a': np.array([0.0, 1.0]), 'sign': '<=', 'b': 0.0})
        # If free, do nothing
    else:
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
        opt_label = "Z (max)" if opt_type.upper() == 'MAX' else "Z (min)"
        opt_val_disp = "-∞" if opt_type.upper() == 'MAX' else "+∞"
        print(f"Hàm mục tiêu {opt_label} = {opt_val_disp}")
        draw_pts = raw_intersections 
    else:
        feasible_pts = sort_clockwise(feasible_pts)
        z_values = np.dot(feasible_pts, c)
        best_z = np.max(z_values) if opt_type == 'MAX' else np.min(z_values)
        optimal_pts = feasible_pts[np.abs(z_values - best_z) < 1e-5]
        is_unbounded = any(abs(abs(pt[i]) - M) < 1e-5 for pt in optimal_pts for i in range(2))
        
        if is_unbounded:
            status = 'UNBOUNDED'
            inf_val = "+∞" if opt_type.upper() == 'MAX' else "-∞"
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
        eqs, basic, non_basic, prob_type, num_vars, num_constraints, geo_constraints, var_signs, z_inputs = parse_input()
        c_np, A_ub, b_ub, A_eq, b_eq = convert_data_for_other_methods(geo_constraints, num_vars, prob_type, z_inputs)
        
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
                solve_dictionary(eqs, basic, non_basic, prob_type, method="Don hinh", var_signs=var_signs, num_vars=num_vars)
            elif choice == '2':
                solve_dictionary(eqs, basic, non_basic, prob_type, method="Bland", var_signs=var_signs, num_vars=num_vars)
            elif choice == '3':
                solve_scipy_2phase(c_np, A_ub, b_ub, A_eq, b_eq, prob_type, var_signs=var_signs)
            elif choice == '4':
                if num_vars != 2:
                    print(f"\n[Lỗi] Phương pháp hình học chỉ áp dụng cho bài toán có 2 ẩn số. Bài toán hiện tại có {num_vars} ẩn.")
                else:
                    solve_geometry(prob_type, c_np, geo_constraints, var_signs=var_signs)
            elif choice == '5':
                print("\n>>> BẮT ĐẦU CHẠY THỬ NGHIỆM ĐỒNG THỜI 4 PHƯƠNG PHÁP <<<")
                solve_dictionary(eqs, basic, non_basic, prob_type, method="Don hinh", var_signs=var_signs, num_vars=num_vars)
                solve_dictionary(eqs, basic, non_basic, prob_type, method="Bland", var_signs=var_signs, num_vars=num_vars)
                solve_scipy_2phase(c_np, A_ub, b_ub, A_eq, b_eq, prob_type, var_signs=var_signs)
                if num_vars == 2:
                    solve_geometry(prob_type, c_np, geo_constraints, var_signs=var_signs)
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