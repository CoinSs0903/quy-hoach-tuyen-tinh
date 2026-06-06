import os
from flask import Flask, request, jsonify, render_template
from solver_web import solve_all_methods

app = Flask(__name__, 
            template_folder='templates',
            static_folder='static',
            static_url_path='/static')

# Ensure directories exist
os.makedirs('templates', exist_ok=True)
os.makedirs('static', exist_ok=True)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/api/solve', methods=['POST'])
def api_solve():
    try:
        data = request.get_json()
        if not data:
            return jsonify({"success": False, "error": "Dữ liệu yêu cầu không hợp lệ (Missing JSON)"}), 400
            
        prob_type = data.get('prob_type', 'min')
        z_coeffs = data.get('z_coeffs', [])
        constraints = data.get('constraints', [])
        
        # Validation
        if not z_coeffs:
            return jsonify({"success": False, "error": "Hệ số hàm mục tiêu Z không được để trống!"}), 400
        if not constraints:
            return jsonify({"success": False, "error": "Vui lòng nhập ít nhất một ràng buộc!"}), 400
            
        # Parse inputs to floats/ints for conversion checks
        # z_coeffs should be list of numbers
        # constraints should be list of dict: {coeffs: list, sign: str, rhs: number}
        parsed_z = [float(x) for x in z_coeffs]
        parsed_constraints = []
        for c in constraints:
            parsed_constraints.append({
                "coeffs": [float(x) for x in c['coeffs']],
                "sign": c['sign'],
                "rhs": float(c['rhs'])
            })
            
        # Run solver
        results = solve_all_methods(prob_type, parsed_z, parsed_constraints)
        
        return jsonify({
            "success": True,
            "results": results,
            "num_vars": len(parsed_z),
            "prob_type": prob_type
        })
        
    except ValueError as ve:
        return jsonify({"success": False, "error": f"Lỗi định dạng số liệu: {str(ve)}"}), 400
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({"success": False, "error": f"Hệ thống gặp lỗi: {str(e)}"}), 500

if __name__ == '__main__':
    # Run locally on port 5000
    app.run(host='127.0.0.1', port=5000, debug=True)
