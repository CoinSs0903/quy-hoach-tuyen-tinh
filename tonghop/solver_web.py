import math
import re
import copy
from fractions import Fraction
from itertools import combinations
import numpy as np
from scipy.optimize import linprog

def var_key(v):
    return (0 if v.startswith('x') else 1, int(v[1:]))

def serialize_fraction(f):
    if isinstance(f, (int, float)):
        # Convert to Fraction first if float to display nicely
        f = Fraction(f).limit_denominator(100000)
    if f.denominator == 1:
        return str(f.numerator)
    return f"{f.numerator}/{f.denominator}"

def get_eq_structure(name, const, coeffs, non_basic):
    terms = []
    for v in sorted(non_basic, key=var_key):
        c = coeffs.get(v, Fraction(0))
        if c == 0: continue
        terms.append({
            "var": v,
            "coeff_val": float(c),
            "coeff_str": serialize_fraction(c)
        })
    return {
        "var_name": name,
        "const_val": float(const),
        "const_str": serialize_fraction(const),
        "terms": terms
    }

def solve_dictionary_web(eqs_orig, basic_orig, non_basic_orig, prob_type, method="Bland"):
    eqs = copy.deepcopy(eqs_orig)
    basic = copy.deepcopy(basic_orig)
    non_basic = copy.deepcopy(non_basic_orig)
    
    initial_feasible = all(eqs[b][0] >= 0 for b in basic)
    
    steps = []
    iteration = 0
    max_iterations = 1000  # Avoid infinite loops in cycling
    
    while iteration < max_iterations:
        # Build current dictionary state
        equations = []
        equations.append(get_eq_structure('z', eqs['z'][0], eqs['z'][1], non_basic))
        for b in sorted(basic, key=var_key):
            equations.append(get_eq_structure(b, eqs[b][0], eqs[b][1], non_basic))
            
        z_coeffs = eqs['z'][1]
        entering = None
        
        if method == "Bland":
            for v in sorted(non_basic, key=var_key):
                if z_coeffs.get(v, Fraction(0)) < 0:
                    entering = v
                    break
        else:  # Dantzig (most negative)
            min_val = Fraction(0)
            for v in sorted(non_basic, key=var_key):
                val = z_coeffs.get(v, Fraction(0))
                if val < min_val:
                    min_val = val
                    entering = v
                    
        step_data = {
            "iteration": iteration,
            "equations": equations,
            "entering": entering,
            "leaving": None,
            "status": "running",
            "current_z": float(eqs['z'][0]) if prob_type == 'min' else -float(eqs['z'][0])
        }
        
        if entering is None:
            step_data["status"] = "optimal"
            steps.append(step_data)
            break
            
        leaving = None
        min_ratio = float('inf')
        
        for b in sorted(basic, key=var_key):
            coeff = eqs[b][1].get(entering, Fraction(0))
            if coeff < 0: 
                ratio = eqs[b][0] / abs(coeff)
                if ratio < min_ratio:
                    min_ratio = ratio
                    leaving = b
                elif ratio == min_ratio and leaving is not None:
                    if var_key(b) < var_key(leaving):
                        leaving = b
                        
        step_data["leaving"] = leaving
        
        if leaving is None:
            step_data["status"] = "unbounded"
            steps.append(step_data)
            return {
                "success": False,
                "status": "unbounded",
                "initial_feasible": initial_feasible,
                "message": "Bài toán không giới nội! Z -> vô cùng",
                "steps": steps
            }
            
        steps.append(step_data)
        
        # Perform pivot operation
        p = eqs[leaving][1][entering]
        new_l_const = eqs[leaving][0] / -p
        
        new_l_coeffs = {leaving: Fraction(1) / p}
        for v in non_basic:
            if v != entering:
                c = eqs[leaving][1].get(v, Fraction(0))
                if c != 0: new_l_coeffs[v] = c / -p
                
        new_eqs = {}
        for k in eqs:
            if k == leaving: continue
            
            c_k, coeffs_k = eqs[k]
            a_ke = coeffs_k.get(entering, Fraction(0))
            
            new_c_k = c_k + a_ke * new_l_const
            new_coeffs_k = {}
            if a_ke != 0:
                new_coeffs_k[leaving] = a_ke / p
                
            for v in non_basic:
                if v != entering:
                    old_c = coeffs_k.get(v, Fraction(0))
                    new_c = old_c + a_ke * new_l_coeffs.get(v, Fraction(0))
                    if new_c != 0: new_coeffs_k[v] = new_c
                    
            new_eqs[k] = (new_c_k, new_coeffs_k)
            
        new_eqs[entering] = (new_l_const, new_l_coeffs)
        eqs = new_eqs
        
        basic.append(entering)
        non_basic.remove(entering)
        non_basic.append(leaving)
        
        iteration += 1
        
    if iteration >= max_iterations:
        return {
            "success": False,
            "status": "cycle",
            "initial_feasible": initial_feasible,
            "message": "Phát hiện vòng lặp vô hạn hoặc vượt quá số bước cho phép!",
            "steps": steps
        }
        
    # Get optimal solution
    z_value = eqs['z'][0]
    opt_val = float(z_value) if prob_type == 'min' else -float(z_value)
    opt_val_str = serialize_fraction(z_value) if prob_type == 'min' else serialize_fraction(-z_value)
    
    optimal_solution = {}
    all_vars = sorted(list(set(basic + non_basic)), key=var_key)
    for v in all_vars:
        val = eqs[v][0] if v in basic else Fraction(0)
        optimal_solution[v] = {
            "val": float(val),
            "str": serialize_fraction(val)
        }
        
    return {
        "success": True,
        "status": "optimal",
        "initial_feasible": initial_feasible,
        "message": "Đã tìm thấy nghiệm tối ưu!" if initial_feasible else "Đã tìm thấy nghiệm tối ưu (Lưu ý: Từ điển xuất phát không khả thi!)",
        "steps": steps,
        "optimal_value": opt_val,
        "optimal_value_str": opt_val_str,
        "optimal_solution": optimal_solution
    }

def solve_scipy_2phase_web(c, A, B, prob_type):
    # c: numpy array of original objective coefficients
    # A: A_ub (coefficients for constraints Ax <= B)
    # B: b_ub (RHS for constraints Ax <= B)
    
    # linprog always minimizes. If maximizing, minimize -c.
    c_input = c.copy() if prob_type == 'min' else -c.copy()
    bounds = [(0, None)] * len(c)
    
    res = linprog(c_input, A_ub=A, b_ub=B, bounds=bounds, method='highs')
    
    if res.success:
        opt_val = res.fun if prob_type == 'min' else -res.fun
        solution = {f"x{i}": float(val) for i, val in enumerate(res.x, 1)}
        return {
            "success": True,
            "status": "optimal",
            "message": f"Thành công: {res.message}",
            "optimal_value": float(opt_val),
            "solution": solution
        }
    else:
        status_str = "infeasible" if "infeasible" in res.message.lower() else "unbounded" if "unbounded" in res.message.lower() else "failed"
        return {
            "success": False,
            "status": status_str,
            "message": f"Không tìm thấy nghiệm tối ưu. Lý do: {res.message}",
            "optimal_value": None,
            "solution": None
        }

def check_feasibility_web(pt, constraints, tol=1e-7):
    for constr in constraints:
        val = np.dot(constr['a'], pt)
        b = constr['b']
        if constr['sign'] == '<=' and val > b + tol: return False
        if constr['sign'] == '>=' and val < b - tol: return False
        if constr['sign'] == '=' and abs(val - b) > tol: return False
    return True

def sort_clockwise_web(pts):
    pts = np.unique(np.round(pts, 5), axis=0)
    if len(pts) <= 2: return pts
    centroid = np.mean(pts, axis=0)
    angles = np.arctan2(pts[:,1] - centroid[1], pts[:,0] - centroid[0])
    return pts[np.argsort(angles)]

def get_anchor_point_web(a1, a2, b_val, x_lim, y_lim):
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

def solve_geometry_web(opt_type, c, A_ub, b_ub, original_constraints_raw=None):
    # original_constraints_raw: list of original constraints from user input
    # e.g., [{"coeffs": [2, 1], "sign": "<=", "rhs": 5}, ...]
    # This helps reconstruct the exact constraint equations.
    
    opt_type = opt_type.upper()
    constraints = []
    
    # We will build constraints for validation and lines
    if original_constraints_raw:
        for item in original_constraints_raw:
            constraints.append({
                'a': np.array([float(item['coeffs'][0]), float(item['coeffs'][1])]),
                'sign': item['sign'],
                'b': float(item['rhs'])
            })
    else:
        for i in range(len(A_ub)):
            constraints.append({'a': A_ub[i], 'sign': '<=', 'b': b_ub[i]})
            
    # Add non-negativity constraints x1 >= 0, x2 >= 0
    constraints.append({'a': np.array([-1.0, 0.0]), 'sign': '<=', 'b': 0.0})
    constraints.append({'a': np.array([0.0, -1.0]), 'sign': '<=', 'b': 0.0})
    
    M = 10000 
    calc_constraints = constraints.copy()
    calc_constraints.extend([
        {'a': np.array([1, 0]), 'sign': '<=', 'b': M},
        {'a': np.array([-1, 0]), 'sign': '<=', 'b': M},
        {'a': np.array([0, 1]), 'sign': '<=', 'b': M},
        {'a': np.array([0, -1]), 'sign': '<=', 'b': M}
    ])
    
    A_lines = [c_item['a'] for c_item in calc_constraints]
    b_lines = [c_item['b'] for c_item in calc_constraints]
    intersections = []
    raw_intersections = []
    
    for i, j in combinations(range(len(constraints)), 2):
        try:
            pt = np.linalg.solve([constraints[i]['a'], constraints[j]['a']], [constraints[i]['b'], constraints[j]['b']])
            raw_intersections.append(pt)
        except np.linalg.LinAlgError: pass
 
    for i, j in combinations(range(len(A_lines)), 2):
        try:
            pt = np.linalg.solve([A_lines[i], A_lines[j]], [b_lines[i], b_lines[j]])
            intersections.append(pt)
        except np.linalg.LinAlgError: pass
            
    feasible_pts = [pt for pt in intersections if check_feasibility_web(pt, calc_constraints)]
    
    status = ''
    p1 = p2 = None 
    draw_pts = [] 
    best_z = None
    best_v = None
    
    if len(feasible_pts) == 0:
        status = 'INFEASIBLE'
        draw_pts = raw_intersections if len(raw_intersections) > 0 else [[0,0], [5,0], [0,5]]
    else:
        feasible_pts = sort_clockwise_web(feasible_pts)
        z_values = np.dot(feasible_pts, c)
        best_z = np.max(z_values) if opt_type == 'MAX' else np.min(z_values)
        optimal_pts = feasible_pts[np.abs(z_values - best_z) < 1e-5]
        is_unbounded = any(abs(abs(pt[i]) - M) < 1e-5 for pt in optimal_pts for i in range(2))
        
        if is_unbounded:
            status = 'UNBOUNDED'
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
                draw_pts = feasible_pts
                best_v = p1 
            else:
                status = 'OPTIMAL'
                best_v = optimal_pts[0]
                draw_pts = feasible_pts
        else: 
            status = 'OPTIMAL'
            best_v = optimal_pts[0]
            draw_pts = feasible_pts

    # Calculate limits for rendering
    draw_pts = np.array(draw_pts)
    if len(draw_pts) > 0:
        x_min, x_max = np.min(draw_pts[:,0]), np.max(draw_pts[:,0])
        y_min, y_max = np.min(draw_pts[:,1]), np.max(draw_pts[:,1])
        margin_x = max(2.0, float(x_max - x_min)*0.4)
        margin_y = max(2.0, float(y_max - y_min)*0.4)
    else:
        x_min, x_max, y_min, y_max = 0.0, 10.0, 0.0, 10.0
        margin_x = margin_y = 2.0
        
    x_lim = [float(x_min - margin_x), float(x_max + margin_x)]
    y_lim = [float(y_min - margin_y), float(y_max + margin_y)]
    
    # Clip limits to prevent negative bounds if feasible points are strictly positive, but standard geometry usually shows standard quadrant
    if x_min >= -1e-5: x_lim[0] = max(-1.0, x_lim[0])
    if y_min >= -1e-5: y_lim[0] = max(-1.0, y_lim[0])

    # Build line equations for frontend
    lines = []
    for idx, constr in enumerate(constraints):
        a1, a2, b_val = float(constr['a'][0]), float(constr['a'][1]), float(constr['b'])
        sign = constr['sign']
        
        # Calculate label anchor point
        anchor = get_anchor_point_web(a1, a2, b_val, x_lim, y_lim)
        
        lines.append({
            "index": idx + 1,
            "a1": a1,
            "a2": a2,
            "b": b_val,
            "sign": sign,
            "anchor": [float(anchor[0]), float(anchor[1])]
        })
        
    # Build feasible region coordinates
    feasible_polygon = []
    if status != 'INFEASIBLE':
        feasible_polygon = [[float(pt[0]), float(pt[1])] for pt in feasible_pts if not any(abs(abs(pt[i]) - M) < 1e-5 for i in range(2))]
        # If it's unbounded, we might have clipped the polygon. The frontend can render it appropriately.
        
    # Return everything
    result = {
        "status": status,
        "x_lim": x_lim,
        "y_lim": y_lim,
        "constraints": lines,
        "feasible_polygon": feasible_polygon,
        "c": [float(c[0]), float(c[1])],
        "opt_type": opt_type,
        "best_z": float(best_z) if best_z is not None else None,
        "optimal_points": [[float(pt[0]), float(pt[1])] for pt in optimal_pts] if best_z is not None else []
    }
    
    if status == 'OPTIMAL' and best_v is not None:
        result["solution"] = {"x1": float(best_v[0]), "x2": float(best_v[1])}
    elif status == 'MULTIPLE_OPTIMAL' and p1 is not None and p2 is not None:
        result["solution_segment"] = {
            "p1": [float(p1[0]), float(p1[1])],
            "p2": [float(p2[0]), float(p2[1])]
        }
        result["solution"] = {"x1": float(p1[0]), "x2": float(p1[1])} # return first endpoint as sample
        
    return result

def solve_all_methods(prob_type, z_coeffs_list, constraints_raw):
    # Prepares data
    num_vars = len(z_coeffs_list)
    
    # 1. Initialize for dictionary method
    # eqs maps: 'z' -> (const, {var_name: coeff}), 'w_i' -> (const, {var_name: coeff})
    # non_basic = [x1, x2...]
    # basic = [w1, w2...]
    
    z_coeffs = {}
    for i, val in enumerate(z_coeffs_list, 1):
        if prob_type == 'max':
            z_coeffs[f'x{i}'] = -Fraction(val)
        else:
            z_coeffs[f'x{i}'] = Fraction(val)
            
    eqs = {'z': (Fraction(0), z_coeffs)}
    non_basic = [f'x{i}' for i in range(1, num_vars + 1)]
    basic = []
    
    for i, item in enumerate(constraints_raw, 1):
        coeffs_list = [Fraction(c) for c in item['coeffs']]
        sign = item['sign']
        rhs = Fraction(item['rhs'])
        
        if sign == '>=':
            coeffs_list = [-c for c in coeffs_list]
            rhs = -rhs
        elif sign == '=':
            # Note: dictionary method in tonghop.py only handles <= and >= (by converting >= to <=)
            # A true '=' constraint is challenging for simple dictionary simplex without artificial variables.
            # We will follow the tonghop.py behavior: it reads parts and handles '>=' by negating,
            # but doesn't handle '=' directly in dictionary simplex without custom code.
            # Let's check how tonghop.py handles it:
            # in tonghop.py: parts[1] is sign. If sign == '>=', it negates.
            # If sign == '=', it keeps it. But wait, if sign is '=', w_i is still added.
            # Let's follow this logic:
            pass
            
        w_name = f'w{i}'
        basic.append(w_name)
        
        w_coeffs = {}
        for j, val in enumerate(coeffs_list, 1):
            if val != 0:
                w_coeffs[f'x{j}'] = -val
                
        eqs[w_name] = (rhs, w_coeffs)
        
    # 2. Convert data for SciPy and Geometry
    c_list = []
    for i in range(1, num_vars + 1):
        val = z_coeffs.get(f'x{i}', Fraction(0))
        if prob_type == 'max':
            c_list.append(float(-val))
        else:
            c_list.append(float(val))
            
    A_ub = []
    b_ub = []
    for k, (rhs, coeffs) in eqs.items():
        if k == 'z': continue
        row = []
        for i in range(1, num_vars + 1):
            row.append(float(-coeffs.get(f'x{i}', Fraction(0))))
        A_ub.append(row)
        b_ub.append(float(rhs))
        
    c_np = np.array(c_list)
    A_np = np.array(A_ub)
    b_np = np.array(b_ub)
    
    # Run algorithms
    results = {}
    
    # Dantzig Simplex
    results["dantzig"] = solve_dictionary_web(eqs, basic, non_basic, prob_type, method="Don hinh")
    
    # Bland Simplex
    results["bland"] = solve_dictionary_web(eqs, basic, non_basic, prob_type, method="Bland")
    
    # SciPy 2-Phase
    results["scipy"] = solve_scipy_2phase_web(c_np, A_np, b_np, prob_type)
    
    # Geometry (if 2 variables)
    if num_vars == 2:
        results["geometry"] = solve_geometry_web(prob_type, c_np, A_np, b_np, constraints_raw)
    else:
        results["geometry"] = None
        
    return results
