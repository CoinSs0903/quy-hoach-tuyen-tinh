import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
import math

def get_user_input():
    print("=== GIẢI QHTT: PHƯƠNG PHÁP TRƯỢT HÀM MỤC TIÊU ===")
    opt_type = input("Loại bài toán (MAX hoặc MIN): ").strip().upper()
    if opt_type not in ['MAX', 'MIN']:
        print("Lỗi: Chỉ chấp nhận MAX hoặc MIN.")
        return None
        
    c_input = input("Nhập ma trận cT (VD: 3 2): ")
    c = np.array([float(val) for val in c_input.split()])
    
    n_constraints = int(input("Nhập tổng số lượng ràng buộc (tính cả x1, x2 >= 0): "))
    
    constraints_info = []
    print("\nNhập ràng buộc theo dạng: a1 a2 dấu b (VD: 2 1 <= 10 hoặc 1 0 >= 0)")
    for i in range(n_constraints):
        constr = input(f"Ràng buộc {i+1}: ").strip().split()
        a1, a2 = float(constr[0]), float(constr[1])
        sign = constr[2]
        b_val = float(constr[3])
        constraints_info.append({'a': np.array([a1, a2]), 'sign': sign, 'b': b_val})
        
    return opt_type, c, constraints_info

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

def main():
    inputs = get_user_input()
    if not inputs: return
    opt_type, c, constraints = inputs
    
    M = 10000 
    calc_constraints = constraints.copy()
    calc_constraints.extend([
        {'a': np.array([1, 0]), 'sign': '<=', 'b': M},
        {'a': np.array([-1, 0]), 'sign': '<=', 'b': M},
        {'a': np.array([0, 1]), 'sign': '<=', 'b': M},
        {'a': np.array([0, -1]), 'sign': '<=', 'b': M}
    ])
    
    A_lines = [c['a'] for c in calc_constraints]
    b_lines = [c['b'] for c in calc_constraints]
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
            
    feasible_pts = [pt for pt in intersections if check_feasibility(pt, calc_constraints)]
    
    print("\n" + "="*55)
    status = ''
    p1 = p2 = None 
    draw_pts = [] 
    
    if len(feasible_pts) == 0:
        status = 'INFEASIBLE'
        print("KẾT LUẬN: BÀI TOÁN VÔ NGHIỆM (Infeasible)")
        print(">> Các ràng buộc mâu thuẫn, không có vùng thỏa mãn chung.")
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
            print("KẾT LUẬN: BÀI TOÁN KHÔNG GIỚI NỘI (Unbounded)")
            print(f">> Miền nghiệm mở ra vô cực. Giá trị tối ưu: {inf_val}")
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
                print("KẾT LUẬN: BÀI TOÁN CÓ VÔ SỐ NGHIỆM")
                print(">> Đường mức hàm mục tiêu trượt trùng khít với 1 cạnh của miền nghiệm.")
                print(f">> Giá trị {opt_type} z* = {best_z:.2f}")
                print(f">> Tập nghiệm là đoạn thẳng nối 2 điểm: A({p1[0]:.2f}, {p1[1]:.2f}) và B({p2[0]:.2f}, {p2[1]:.2f})")
                
                print("\n[Biểu diễn phương trình nghiệm]")
                if abs(p2[0] - p1[0]) < 1e-5:
                    print(f"  x1 = {p1[0]:.2f}")
                    print(f"  x2 thuộc [{min(p1[1], p2[1]):.2f}, {max(p1[1], p2[1]):.2f}]")
                elif abs(p2[1] - p1[1]) < 1e-5:
                    print(f"  x2 = {p1[1]:.2f}")
                    print(f"  x1 thuộc [{min(p1[0], p2[0]):.2f}, {max(p1[0], p2[0]):.2f}]")
                else:
                    slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
                    intercept = p1[1] - slope * p1[0]
                    sign_i = "+" if intercept >= 0 else "-"
                    print(f"  x2 = {slope:.2f} * x1  {sign_i}  {abs(intercept):.2f}")
                    print(f"  x1 thuộc [{min(p1[0], p2[0]):.2f}, {max(p1[0], p2[0]):.2f}]")
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
            print(f">> Giao điểm tối ưu: x1 = {best_v[0]:.2f}, x2 = {best_v[1]:.2f}")
            print(f">> Giá trị hàm mục tiêu {opt_type} z = {best_z:.2f}")
            draw_pts = feasible_pts
            
    print("="*55)

    # VẼ ĐỒ THỊ 
    fig, ax = plt.subplots(figsize=(10, 8))
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
        a1, a2, b_val, sign = constr['a'][0], constr['a'][1], constr['b'], constr['sign']
        if abs(a2) > 1e-7:
            x_vals = np.array(x_lim_draw)
            ax.plot(x_vals, (b_val - a1 * x_vals) / a2, color='black', lw=1.5, alpha=0.6)
        else:
            ax.axvline(x=b_val/a1, color='black', lw=1.5, alpha=0.6)

        nx, ny = a1, a2
        if sign == '<=': nx, ny = -nx, -ny
        elif sign == '=': nx, ny = 0, 0
        norm = math.hypot(nx, ny)
        if norm > 1e-7: nx, ny = (nx / norm) * arrow_len, (ny / norm) * arrow_len
            
        anchor = get_anchor_point(a1, a2, b_val, x_lim_draw, y_lim_draw)
        ax.text(anchor[0], anchor[1], f'({idx+1})', color='black', fontsize=11, fontweight='bold', 
                ha='center', va='center', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', boxstyle='circle'))
        
        if norm > 1e-7:
            ax.annotate('', xy=(anchor[0] + nx, anchor[1] + ny), xytext=(anchor[0], anchor[1]),
                        arrowprops=dict(arrowstyle="->,head_width=0.4,head_length=0.6", color='blue', lw=1.5))

    if status != 'INFEASIBLE' and len(feasible_pts) > 2:
        polygon = plt.Polygon(feasible_pts, closed=True, color='cyan', alpha=0.4, label='Miền nghiệm')
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
            
        draw_objective_line(ax, c, best_z, x_lim_draw, color='red', linestyle='--', lw=2.5, label=f'Đường mức z*={best_z:.2f}')
        
        if status == 'MULTIPLE_OPTIMAL':
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color='magenta', lw=5, zorder=6, label='Tập vô số nghiệm')
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], marker='*', color='red', markersize=15, linestyle='None', zorder=7)
            vector_origin = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2)
        else:
            ax.plot(best_v[0], best_v[1], marker='*', color='red', markersize=15, zorder=5)
            vector_origin = (best_v[0], best_v[1])
        
        ax.annotate('', xy=(vector_origin[0] + zx, vector_origin[1] + zy), xytext=(vector_origin[0], vector_origin[1]),
                    arrowprops=dict(arrowstyle="simple,head_width=0.6,head_length=0.8", color='red', lw=1.5))
        ax.text(vector_origin[0] + zx, vector_origin[1] + zy, ' Hướng z', color='red', fontsize=12, fontweight='bold', va='bottom')

    elif status == 'UNBOUNDED':
        center_x, center_y = np.mean(draw_pts[:,0]), np.mean(draw_pts[:,1])
        ax.annotate('', xy=(center_x + zx*2, center_y + zy*2), xytext=(center_x, center_y),
                    arrowprops=dict(arrowstyle="simple,head_width=0.8,head_length=1.0", color='red', lw=2.0))
        ax.text(center_x + zx*2, center_y + zy*2, f' Tiến về {inf_val}', color='red', fontsize=12, fontweight='bold', va='bottom')

    elif status == 'INFEASIBLE':
        ax.text(sum(x_lim_draw)/2, sum(y_lim_draw)/2, 'VÔ NGHIỆM\n(Không có miền chung)', 
                color='red', fontsize=16, fontweight='bold', ha='center', va='center', alpha=0.5,
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='red'))

    ax.set_xlim(x_lim_draw)
    ax.set_ylim(y_lim_draw)
    ax.axhline(0, color='black', lw=1)
    ax.axvline(0, color='black', lw=1)
    ax.set_xlabel(r'$x_1$')
    ax.set_ylabel(r'$x_2$')
    
    title = f'Phương Pháp Trượt: {opt_type} z = {c[0]}x1 + {c[1]}x2'
    if status != 'OPTIMAL': title += f" ({status})"
    ax.set_title(title, fontweight='bold', fontsize=14)
    
    if status in ['OPTIMAL', 'MULTIPLE_OPTIMAL']: ax.legend(loc='upper right')
    ax.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()