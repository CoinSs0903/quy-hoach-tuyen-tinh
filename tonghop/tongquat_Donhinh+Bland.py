from fractions import Fraction
import re

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
    print("=== NHẬP DỮ LIỆU BÀI TOÁN TỔNG QUÁT ===")
    
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
            
    
    num_constraints = int(input("\nNhập số lượng ràng buộc: "))
    
    eqs = {'z': (Fraction(0), z_coeffs)}
    non_basic = [f'x{i}' for i in range(1, num_vars + 1)]
    basic = []
    
    print("\nNhập các ràng buộc (Ví dụ: 2x1 + 2x2 <= 9 nhập là: 2 2 <= 9)")
    for i in range(1, num_constraints + 1):
        line = input(f"Ràng buộc {i}: ").strip()
        
        parts = re.split(r'\s*(<=|>=|=)\s*', line)
        
        coeffs_list = list(map(Fraction, parts[0].split()))
        sign = parts[1]
        rhs = Fraction(parts[2])
        
       
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
        
    return eqs, basic, non_basic, prob_type

def solve_dictionary(eqs_orig, basic_orig, non_basic_orig, prob_type, method="Bland"):
    import copy
    eqs = copy.deepcopy(eqs_orig)
    basic = copy.deepcopy(basic_orig)
    non_basic = copy.deepcopy(non_basic_orig)
    
    print(f"\n{'='*60}")
    print(f"GIẢI BẰNG PHƯƠNG PHÁP: {method.upper()} ({prob_type.upper()})")
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
            print("\n>> Bài toán không giới nội (Vô cùng)! Z -> +vô cùng")
            return
            
        print(f"\n=> Xoay: Biến VÀO = {entering}, Biến RA = {leaving}")
        

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
        
        basic.append(entering)
        non_basic.remove(entering)
        non_basic.append(leaving)
        
        iteration += 1

    print(f"\n{'*'*20} KẾT QUẢ TỐI ƯU {'*'*20}")
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


if __name__ == "__main__":
    try:
        # 1. Nhập bài toán tổng quát từ bàn phím
        eqs, basic, non_basic, prob_type = parse_input()
        
        # 2. Giải bằng cả 2 phương pháp để so sánh
        solve_dictionary(eqs, basic, non_basic, prob_type, method="Bland")
        solve_dictionary(eqs, basic, non_basic, prob_type, method="Don hinh")
    except Exception as e:
        print(f"\n[Lỗi] Vui lòng kiểm tra lại định dạng nhập vào! Chi tiết: {e}")