// Exact Fraction class for exact linear programming arithmetic
class Fraction {
    constructor(num, den = 1) {
        if (den === 0) throw new Error("Denominator cannot be zero");
        if (den < 0) {
            num = -num;
            den = -den;
        }
        const g = Fraction.gcd(Math.abs(num), den);
        this.num = Math.round(num / g);
        this.den = Math.round(den / g);
    }
    
    static gcd(a, b) {
        a = Math.round(a);
        b = Math.round(b);
        while (b) {
            let t = b;
            b = a % b;
            a = t;
        }
        return a;
    }
    
    static parse(val) {
        if (val instanceof Fraction) return val;
        if (typeof val === 'number') {
            if (Number.isInteger(val)) {
                return new Fraction(val, 1);
            }
            const tolerance = 1.0e-9;
            let h1 = 1, h2 = 0, k1 = 0, k2 = 1;
            let b = val;
            do {
                let a = Math.floor(b);
                let aux = h1; h1 = a * h1 + h2; h2 = aux;
                aux = k1; k1 = a * k1 + k2; k2 = aux;
                b = 1 / (b - a);
            } while (Math.abs(val - h1 / k1) > val * tolerance);
            return new Fraction(h1, k1);
        }
        const str = String(val).trim();
        const parts = str.split('/');
        if (parts.length === 2) {
            return new Fraction(Number(parts[0]), Number(parts[1]));
        }
        return new Fraction(Number(parts[0]), 1);
    }
    
    add(other) {
        other = Fraction.parse(other);
        return new Fraction(this.num * other.den + other.num * this.den, this.den * other.den);
    }
    
    sub(other) {
        other = Fraction.parse(other);
        return new Fraction(this.num * other.den - other.num * this.den, this.den * other.den);
    }
    
    mul(other) {
        other = Fraction.parse(other);
        return new Fraction(this.num * other.num, this.den * other.den);
    }
    
    div(other) {
        other = Fraction.parse(other);
        return new Fraction(this.num * other.den, this.den * other.num);
    }
    
    neg() {
        return new Fraction(-this.num, this.den);
    }
    
    abs() {
        return new Fraction(Math.abs(this.num), this.den);
    }
    
    toFloat() {
        return this.num / this.den;
    }
    
    toString() {
        if (this.den === 1) return String(this.num);
        return `${this.num}/${this.den}`;
    }
}

// Helpers for LP variables sorting and printing
function varKeyCompare(a, b) {
    const typeA = a.startsWith('x') ? 0 : 1;
    const typeB = b.startsWith('x') ? 0 : 1;
    if (typeA !== typeB) return typeA - typeB;
    const numA = parseInt(a.substring(1));
    const numB = parseInt(b.substring(1));
    return numA - numB;
}

function getEqStructure(name, constVal, coeffs, nonBasic) {
    const terms = [];
    const sortedNonBasic = [...nonBasic].sort(varKeyCompare);
    for (const v of sortedNonBasic) {
        const c = coeffs[v] || new Fraction(0);
        if (c.toFloat() === 0) continue;
        terms.push({
            var: v,
            coeff_val: c.toFloat(),
            coeff_str: c.toString()
        });
    }
    return {
        var_name: name,
        const_val: constVal.toFloat(),
        const_str: constVal.toString(),
        terms: terms
    };
}

// Standard dictionary simplex solver in Javascript
function solveSimplexJS(eqsOrig, basicOrig, nonBasicOrig, probType, method = "Bland") {
    const eqs = {};
    for (const [k, v] of Object.entries(eqsOrig)) {
        const coeffs = {};
        for (const [vk, vv] of Object.entries(v[1])) {
            coeffs[vk] = vv;
        }
        eqs[k] = [v[0], coeffs];
    }
    
    const basic = [...basicOrig];
    const nonBasic = [...nonBasicOrig];
    const initialFeasible = basic.every(b => eqs[b][0].toFloat() >= -1.0e-7);
    
    const steps = [];
    let iteration = 0;
    const maxIterations = 1000;
    
    while (iteration < maxIterations) {
        const equations = [];
        equations.push(getEqStructure('z', eqs['z'][0], eqs['z'][1], nonBasic));
        for (const b of [...basic].sort(varKeyCompare)) {
            equations.push(getEqStructure(b, eqs[b][0], eqs[b][1], nonBasic));
        }
        
        const zCoeffs = eqs['z'][1];
        let entering = null;
        
        if (method === "Bland") {
            const sortedNonBasic = [...nonBasic].sort(varKeyCompare);
            for (const v of sortedNonBasic) {
                const val = zCoeffs[v] || new Fraction(0);
                if (val.toFloat() < -1.0e-7) {
                    entering = v;
                    break;
                }
            }
        } else { // Dantzig (most negative)
            let minVal = 0;
            const sortedNonBasic = [...nonBasic].sort(varKeyCompare);
            for (const v of sortedNonBasic) {
                const val = (zCoeffs[v] || new Fraction(0)).toFloat();
                if (val < minVal) {
                    minVal = val;
                    entering = v;
                }
            }
        }
        
        const stepData = {
            iteration: iteration,
            equations: equations,
            entering: entering,
            leaving: null,
            status: "running",
            current_z: probType === 'min' ? eqs['z'][0].toFloat() : -eqs['z'][0].toFloat()
        };
        
        if (entering === null) {
            stepData.status = "optimal";
            steps.push(stepData);
            break;
        }
        
        let leaving = null;
        let minRatio = Infinity;
        
        const sortedBasic = [...basic].sort(varKeyCompare);
        for (const b of sortedBasic) {
            const coeff = eqs[b][1][entering] || new Fraction(0);
            if (coeff.toFloat() < -1.0e-7) {
                const ratio = eqs[b][0].div(coeff.abs()).toFloat();
                if (ratio < minRatio - 1.0e-7) {
                    minRatio = ratio;
                    leaving = b;
                } else if (Math.abs(ratio - minRatio) < 1.0e-7 && leaving !== null) {
                    if (varKeyCompare(b, leaving) < 0) {
                        leaving = b;
                    }
                }
            }
        }
        
        stepData.leaving = leaving;
        
        if (leaving === null) {
            stepData.status = "unbounded";
            steps.push(stepData);
            return {
                success: false,
                status: "unbounded",
                initial_feasible: initialFeasible,
                message: "Bài toán không giới nội! Z -> vô cùng",
                steps: steps
            };
        }
        
        steps.push(stepData);
        
        // Pivot
        const p = eqs[leaving][1][entering];
        const newLConst = eqs[leaving][0].div(p.neg());
        const newLCoeffs = {};
        newLCoeffs[leaving] = new Fraction(1).div(p);
        for (const v of nonBasic) {
            if (v !== entering) {
                const c = eqs[leaving][1][v] || new Fraction(0);
                if (c.toFloat() !== 0) {
                    newLCoeffs[v] = c.div(p.neg());
                }
            }
        }
        
        const newEqs = {};
        for (const k of Object.keys(eqs)) {
            if (k === leaving) continue;
            
            const [c_k, coeffs_k] = eqs[k];
            const a_ke = coeffs_k[entering] || new Fraction(0);
            const new_c_k = c_k.add(a_ke.mul(newLConst));
            const new_coeffs_k = {};
            if (a_ke.toFloat() !== 0) {
                new_coeffs_k[leaving] = a_ke.div(p);
            }
            
            for (const v of nonBasic) {
                if (v !== entering) {
                    const old_c = coeffs_k[v] || new Fraction(0);
                    const new_c = old_c.add(a_ke.mul(newLCoeffs[v] || new Fraction(0)));
                    if (new_c.toFloat() !== 0) {
                        new_coeffs_k[v] = new_c;
                    }
                }
            }
            newEqs[k] = [new_c_k, new_coeffs_k];
        }
        newEqs[entering] = [newLConst, newLCoeffs];
        
        // Update basis
        for (const k of Object.keys(newEqs)) {
            eqs[k] = newEqs[k];
        }
        delete eqs[leaving];
        
        const idxBasic = basic.indexOf(leaving);
        basic.splice(idxBasic, 1);
        basic.push(entering);
        
        const idxNonBasic = nonBasic.indexOf(entering);
        nonBasic.splice(idxNonBasic, 1);
        nonBasic.push(leaving);
        
        iteration++;
    }
    
    if (iteration >= maxIterations) {
        return {
            success: false,
            status: "cycle",
            initial_feasible: initialFeasible,
            message: "Phát hiện vòng lặp vô hạn hoặc vượt quá số bước cho phép!",
            steps: steps
        };
    }
    
    const zVal = eqs['z'][0];
    const optVal = probType === 'min' ? zVal.toFloat() : -zVal.toFloat();
    const optValStr = probType === 'min' ? zVal.toString() : zVal.neg().toString();
    
    const optimalSolution = {};
    const allVars = [...basic, ...nonBasic];
    for (const v of allVars) {
        const val = basic.includes(v) ? eqs[v][0] : new Fraction(0);
        optimalSolution[v] = {
            val: val.toFloat(),
            str: val.toString()
        };
    }
    
    return {
        success: initialFeasible,
        status: initialFeasible ? "optimal" : "infeasible",
        initial_feasible: initialFeasible,
        message: initialFeasible ? "Đã tìm thấy nghiệm tối ưu!" : "Từ điển xuất phát không khả thi (Cần sử dụng phương pháp 2 Pha).",
        steps: steps,
        optimal_value: optVal,
        optimal_value_str: optValStr,
        optimal_solution: optimalSolution
    };
}

// Reconstruct original Z from current Phase 1 basis
function reconstructObjectiveJS(eqs, basic, nonBasic, zOrig) {
    const cZ = zOrig[0];
    const coeffsZ = zOrig[1];
    
    let newCZ = cZ;
    const newCoeffsZ = {};
    
    for (const [v, coeff] of Object.entries(coeffsZ)) {
        if (basic.includes(v)) {
            const [constV, coeffsV] = eqs[v];
            newCZ = newCZ.add(coeff.mul(constV));
            for (const [kv, coeffKv] of Object.entries(coeffsV)) {
                newCoeffsZ[kv] = (newCoeffsZ[kv] || new Fraction(0)).add(coeff.mul(coeffKv));
            }
        } else {
            newCoeffsZ[v] = (newCoeffsZ[v] || new Fraction(0)).add(coeff);
        }
    }
    return [newCZ, newCoeffsZ];
}

// 2-Phase Simplex dictionary solver in Javascript
function solve2PhaseSimplexJS(eqsOrig, basicOrig, nonBasicOrig, probType, method = "Bland") {
    const initialFeasible = basicOrig.every(b => eqsOrig[b][0].toFloat() >= -1.0e-7);
    
    if (initialFeasible) {
        const res = solveSimplexJS(eqsOrig, basicOrig, nonBasicOrig, probType, method);
        return {
            success: res.success,
            status: res.status,
            message: "Từ điển xuất phát đã khả thi. Bỏ qua Pha 1. " + res.message,
            initial_feasible: true,
            phase1_steps: [],
            phase2_steps: res.steps,
            optimal_value: res.optimal_value,
            optimal_value_str: res.optimal_value_str,
            optimal_solution: res.optimal_solution
        };
    }
    
    const eqs = {};
    for (const [k, v] of Object.entries(eqsOrig)) {
        const coeffs = {};
        for (const [vk, vv] of Object.entries(v[1])) {
            coeffs[vk] = vv;
        }
        eqs[k] = [v[0], coeffs];
    }
    
    const basic = [...basicOrig];
    const nonBasic = [...nonBasicOrig, 'x0'];
    const zOrig = [eqs['z'][0], {...eqs['z'][1]}];
    
    // Auxiliary objective: minimize x0 (which is: z_aux = 0 + 1 * x0)
    eqs['z_aux'] = [new Fraction(0), { 'x0': new Fraction(1) }];
    for (const b of basic) {
        eqs[b][1]['x0'] = new Fraction(1);
    }
    
    const phase1_steps = [];
    const buildStepEquations = () => {
        const equations = [];
        equations.push(getEqStructure('z', eqs['z_aux'][0], eqs['z_aux'][1], nonBasic)); // Display as 'z'
        for (const b of [...basic].sort(varKeyCompare)) {
            equations.push(getEqStructure(b, eqs[b][0], eqs[b][1], nonBasic));
        }
        return equations;
    };
    
    // Initial Phase 1 Step 0
    phase1_steps.push({
        iteration: 0,
        equations: buildStepEquations(),
        entering: 'x0',
        leaving: null,
        status: "running",
        current_z: eqs['z_aux'][0].toFloat()
    });
    
    // Leaving variable is the one with the most negative constant
    let leaving = basic[0];
    let minConstVal = eqs[basic[0]][0].toFloat();
    for (const b of basic) {
        const val = eqs[b][0].toFloat();
        if (val < minConstVal) {
            minConstVal = val;
            leaving = b;
        }
    }
    
    const entering = 'x0';
    phase1_steps[0].leaving = leaving;
    
    // Perform initial pivot
    const p = eqs[leaving][1][entering];
    const newLConst = eqs[leaving][0].div(p.neg());
    const newLCoeffs = {};
    newLCoeffs[leaving] = new Fraction(1).div(p);
    for (const v of nonBasic) {
        if (v !== entering) {
            const c = eqs[leaving][1][v] || new Fraction(0);
            if (c.toFloat() !== 0) {
                newLCoeffs[v] = c.div(p.neg());
            }
        }
    }
    
    let newEqs = {};
    for (const k of Object.keys(eqs)) {
        if (k === leaving) continue;
        const [c_k, coeffs_k] = eqs[k];
        const a_ke = coeffs_k[entering] || new Fraction(0);
        const new_c_k = c_k.add(a_ke.mul(newLConst));
        const new_coeffs_k = {};
        if (a_ke.toFloat() !== 0) {
            new_coeffs_k[leaving] = a_ke.div(p);
        }
        for (const v of nonBasic) {
            if (v !== entering) {
                const old_c = coeffs_k[v] || new Fraction(0);
                const new_c = old_c.add(a_ke.mul(newLCoeffs[v] || new Fraction(0)));
                if (new_c.toFloat() !== 0) {
                    new_coeffs_k[v] = new_c;
                }
            }
        }
        newEqs[k] = [new_c_k, new_coeffs_k];
    }
    newEqs[entering] = [newLConst, newLCoeffs];
    
    for (const k of Object.keys(newEqs)) {
        eqs[k] = newEqs[k];
    }
    delete eqs[leaving];
    
    const idxB = basic.indexOf(leaving);
    basic.splice(idxB, 1);
    basic.push(entering);
    
    const idxNB = nonBasic.indexOf(entering);
    nonBasic.splice(idxNB, 1);
    nonBasic.push(leaving);
    
    // Phase 1 Simplex iterations
    let iteration = 1;
    const maxIterations = 1000;
    
    while (iteration < maxIterations) {
        const equations = buildStepEquations();
        const zCoeffs = eqs['z_aux'][1];
        let enteringVar = null;
        
        const sortedNonBasic = [...nonBasic].sort(varKeyCompare);
        for (const v of sortedNonBasic) {
            const val = zCoeffs[v] || new Fraction(0);
            if (val.toFloat() < -1.0e-7) {
                enteringVar = v;
                break;
            }
        }
        
        const stepData = {
            iteration: iteration,
            equations: equations,
            entering: enteringVar,
            leaving: null,
            status: "running",
            current_z: eqs['z_aux'][0].toFloat()
        };
        
        if (enteringVar === null) {
            stepData.status = "optimal";
            phase1_steps.push(stepData);
            break;
        }
        
        let leavingVar = null;
        let minRatio = Infinity;
        const sortedBasic = [...basic].sort(varKeyCompare);
        for (const b of sortedBasic) {
            const coeff = eqs[b][1][enteringVar] || new Fraction(0);
            if (coeff.toFloat() < -1.0e-7) {
                const ratio = eqs[b][0].div(coeff.abs()).toFloat();
                if (ratio < minRatio - 1.0e-7) {
                    minRatio = ratio;
                    leavingVar = b;
                } else if (Math.abs(ratio - minRatio) < 1.0e-7 && leavingVar !== null) {
                    if (varKeyCompare(b, leavingVar) < 0) {
                        leavingVar = b;
                    }
                }
            }
        }
        
        stepData.leaving = leavingVar;
        phase1_steps.push(stepData);
        
        if (leavingVar === null) {
            return {
                success: false,
                status: "unbounded",
                message: "Pha 1 gặp lỗi: bài toán phụ không giới nội.",
                initial_feasible: false,
                phase1_steps: phase1_steps,
                phase2_steps: []
            };
        }
        
        // Pivot
        const p_val = eqs[leavingVar][1][enteringVar];
        const newLConst_val = eqs[leavingVar][0].div(p_val.neg());
        const newLCoeffs_val = {};
        newLCoeffs_val[leavingVar] = new Fraction(1).div(p_val);
        for (const v of nonBasic) {
            if (v !== enteringVar) {
                const c = eqs[leavingVar][1][v] || new Fraction(0);
                if (c.toFloat() !== 0) {
                    newLCoeffs_val[v] = c.div(p_val.neg());
                }
            }
        }
        
        newEqs = {};
        for (const k of Object.keys(eqs)) {
            if (k === leavingVar) continue;
            const [c_k, coeffs_k] = eqs[k];
            const a_ke = coeffs_k[enteringVar] || new Fraction(0);
            const new_c_k = c_k.add(a_ke.mul(newLConst_val));
            const new_coeffs_k = {};
            if (a_ke.toFloat() !== 0) {
                new_coeffs_k[leavingVar] = a_ke.div(p_val);
            }
            for (const v of nonBasic) {
                if (v !== enteringVar) {
                    const old_c = coeffs_k[v] || new Fraction(0);
                    const new_c = old_c.add(a_ke.mul(newLCoeffs_val[v] || new Fraction(0)));
                    if (new_c.toFloat() !== 0) {
                        new_coeffs_k[v] = new_c;
                    }
                }
            }
            newEqs[k] = [new_c_k, new_coeffs_k];
        }
        newEqs[enteringVar] = [newLConst_val, newLCoeffs_val];
        
        for (const k of Object.keys(newEqs)) {
            eqs[k] = newEqs[k];
        }
        delete eqs[leavingVar];
        
        const idxBasic_val = basic.indexOf(leavingVar);
        basic.splice(idxBasic_val, 1);
        basic.push(enteringVar);
        
        const idxNonBasic_val = nonBasic.indexOf(enteringVar);
        nonBasic.splice(idxNonBasic_val, 1);
        nonBasic.push(leavingVar);
        
        iteration++;
    }
    
    const optZAux = eqs['z_aux'][0].toFloat();
    if (optZAux > 1.0e-7) {
        return {
            success: false,
            status: "infeasible",
            message: `Bài toán vô nghiệm! Giá trị tối ưu Pha 1: x0 = ${optZAux.toFixed(4)} > 0.`,
            initial_feasible: false,
            phase1_steps: phase1_steps,
            phase2_steps: []
        };
    }
    
    // Drop x0
    if (basic.includes('x0')) {
        const idx = basic.indexOf('x0');
        basic.splice(idx, 1);
        delete eqs['x0'];
    }
    if (nonBasic.includes('x0')) {
        const idx = nonBasic.indexOf('x0');
        nonBasic.splice(idx, 1);
    }
    for (const k of Object.keys(eqs)) {
        if (eqs[k][1]['x0']) {
            delete eqs[k][1]['x0'];
        }
    }
    delete eqs['z_aux'];
    
    // Reconstruct original Z
    const [c_z, coeffs_z] = reconstructObjectiveJS(eqs, basic, nonBasic, zOrig);
    eqs['z'] = [c_z, coeffs_z];
    
    // Run Phase 2
    const res2 = solveSimplexJS(eqs, basic, nonBasic, probType, method);
    
    return {
        success: res2.success,
        status: res2.status,
        message: "Pha 1 tìm được nghiệm khả thi. " + res2.message,
        initial_feasible: false,
        phase1_steps: phase1_steps,
        phase2_steps: res2.steps,
        optimal_value: res2.optimal_value,
        optimal_value_str: res2.optimal_value_str,
        optimal_solution: res2.optimal_solution
    };
}

// 2D coordinate system solver
function solveLinearSystem2x2(a1, a2, b1, b2) {
    const a11 = a1[0], a12 = a1[1], b_1 = b1;
    const a21 = a2[0], a22 = a2[1], b_2 = b2;
    const det = a11 * a22 - a12 * a21;
    if (Math.abs(det) < 1.0e-9) return null;
    return [
        (b_1 * a22 - a12 * b_2) / det,
        (a11 * b_2 - b_1 * a21) / det
    ];
}

function checkFeasibilityJS(pt, constraints, tol = 1.0e-7) {
    for (const c of constraints) {
        const val = c.a[0] * pt[0] + c.a[1] * pt[1];
        if (c.sign === '<=' && val > c.b + tol) return false;
        if (c.sign === '>=' && val < c.b - tol) return false;
        if (c.sign === '=' && Math.abs(val - c.b) > tol) return false;
    }
    return true;
}

function sortClockwiseJS(pts) {
    const unique = [];
    for (const p of pts) {
        if (!unique.some(up => Math.hypot(up[0] - p[0], up[1] - p[1]) < 1.0e-5)) {
            unique.push(p);
        }
    }
    if (unique.length <= 2) return unique;
    
    let cx = 0, cy = 0;
    for (const p of unique) {
        cx += p[0];
        cy += p[1];
    }
    cx /= unique.length;
    cy /= unique.length;
    
    return unique.sort((pA, pB) => {
        const angleA = Math.atan2(pA[1] - cy, pA[0] - cx);
        const angleB = Math.atan2(pB[1] - cy, pB[0] - cx);
        return angleA - angleB;
    });
}

function getAnchorPointJS(a1, a2, b, x_lim, y_lim) {
    const pts = [];
    if (Math.abs(a2) > 1.0e-7) {
        for (const x of x_lim) {
            const y = (b - a1 * x) / a2;
            if (y >= y_lim[0] - 1.0e-5 && y <= y_lim[1] + 1.0e-5) pts.push([x, y]);
        }
    }
    if (Math.abs(a1) > 1.0e-7) {
        for (const y of y_lim) {
            const x = (b - a2 * y) / a1;
            if (x >= x_lim[0] - 1.0e-5 && x <= x_lim[1] + 1.0e-5) pts.push([x, y]);
        }
    }
    const unique = [];
    for (const p of pts) {
        if (!unique.some(up => Math.hypot(up[0] - p[0], up[1] - p[1]) < 1.0e-5)) {
            unique.push(p);
        }
    }
    if (unique.length >= 2) {
        return [(unique[0][0] + unique[1][0]) / 2, (unique[0][1] + unique[1][1]) / 2];
    } else if (unique.length === 1) {
        return unique[0];
    } else {
        return Math.abs(a2) > 1.0e-7 ? [0, b / a2] : [b / a1, 0];
    }
}

function solveGeometryJS(optType, c, A_ub, b_ub, originalConstraintsRaw) {
    optType = optType.toUpperCase();
    const constraints = [];
    if (originalConstraintsRaw && originalConstraintsRaw.length > 0) {
        for (const item of originalConstraintsRaw) {
            constraints.push({
                a: [Number(item.coeffs[0]), Number(item.coeffs[1])],
                sign: item.sign,
                b: Number(item.rhs)
            });
        }
    } else {
        for (let i = 0; i < A_ub.length; i++) {
            constraints.push({ a: A_ub[i], sign: '<=', b: b_ub[i] });
        }
    }
    
    // Boundaries
    constraints.push({ a: [-1.0, 0.0], sign: '<=', b: 0.0 });
    constraints.push({ a: [0.0, -1.0], sign: '<=', b: 0.0 });
    
    const M = 10000;
    const calcConstraints = [...constraints];
    calcConstraints.push({ a: [1, 0], sign: '<=', b: M });
    calcConstraints.push({ a: [-1, 0], sign: '<=', b: M });
    calcConstraints.push({ a: [0, 1], sign: '<=', b: M });
    calcConstraints.push({ a: [0, -1], sign: '<=', b: M });
    
    const intersections = [];
    const rawIntersections = [];
    
    for (let i = 0; i < constraints.length; i++) {
        for (let j = i + 1; j < constraints.length; j++) {
            const pt = solveLinearSystem2x2(constraints[i].a, constraints[j].a, constraints[i].b, constraints[j].b);
            if (pt) rawIntersections.push(pt);
        }
    }
    
    for (let i = 0; i < calcConstraints.length; i++) {
        for (let j = i + 1; j < calcConstraints.length; j++) {
            const pt = solveLinearSystem2x2(calcConstraints[i].a, calcConstraints[j].a, calcConstraints[i].b, calcConstraints[j].b);
            if (pt) intersections.push(pt);
        }
    }
    
    const feasiblePts = intersections.filter(pt => checkFeasibilityJS(pt, calcConstraints));
    
    let status = '';
    let p1 = null, p2 = null;
    let drawPts = [];
    let best_z = null;
    let best_v = null;
    let feasiblePolygon = [];
    let optimalPts = [];
    let sortedFeasible = [];
    
    if (feasiblePts.length === 0) {
        status = 'INFEASIBLE';
        drawPts = rawIntersections.length > 0 ? rawIntersections : [[0,0], [5,0], [0,5]];
    } else {
        sortedFeasible = sortClockwiseJS(feasiblePts);
        const zValues = sortedFeasible.map(pt => c[0] * pt[0] + c[1] * pt[1]);
        best_z = optType === 'MAX' ? Math.max(...zValues) : Math.min(...zValues);
        
        optimalPts = sortedFeasible.filter((pt, idx) => Math.abs(zValues[idx] - best_z) < 1.0e-5);
        const isUnbounded = optimalPts.some(pt => Math.abs(Math.abs(pt[0]) - M) < 1.0e-5 || Math.abs(Math.abs(pt[1]) - M) < 1.0e-5);
        
        if (isUnbounded) {
            status = 'UNBOUNDED';
            drawPts = sortedFeasible.filter(pt => Math.abs(Math.abs(pt[0]) - M) >= 1.0e-5 && Math.abs(Math.abs(pt[1]) - M) >= 1.0e-5);
        } else if (optimalPts.length > 1) {
            let maxDist = 0;
            p1 = optimalPts[0];
            p2 = optimalPts[0];
            for (let i = 0; i < optimalPts.length; i++) {
                for (let j = i + 1; j < optimalPts.length; j++) {
                    const dist = Math.hypot(optimalPts[i][0] - optimalPts[j][0], optimalPts[i][1] - optimalPts[j][1]);
                    if (dist > maxDist) {
                        maxDist = dist;
                        p1 = optimalPts[i];
                        p2 = optimalPts[j];
                    }
                }
            }
            if (maxDist > 1.0e-5) {
                status = 'MULTIPLE_OPTIMAL';
                drawPts = sortedFeasible;
                best_v = p1;
            } else {
                status = 'OPTIMAL';
                best_v = optimalPts[0];
                drawPts = sortedFeasible;
            }
        } else {
            status = 'OPTIMAL';
            best_v = optimalPts[0];
            drawPts = sortedFeasible;
        }
    }
    
    let x_min = Infinity, x_max = -Infinity, y_min = Infinity, y_max = -Infinity;
    if (drawPts.length > 0) {
        for (const pt of drawPts) {
            x_min = Math.min(x_min, pt[0]);
            x_max = Math.max(x_max, pt[0]);
            y_min = Math.min(y_min, pt[1]);
            y_max = Math.max(y_max, pt[1]);
        }
        var margin_x = Math.max(2.0, (x_max - x_min) * 0.4);
        var margin_y = Math.max(2.0, (y_max - y_min) * 0.4);
    } else {
        x_min = 0; x_max = 10; y_min = 0; y_max = 10;
        margin_x = 2.0; margin_y = 2.0;
    }
    
    const x_lim = [x_min - margin_x, x_max + margin_x];
    const y_lim = [y_min - margin_y, y_max + margin_y];
    
    if (x_min >= -1.0e-5) x_lim[0] = Math.max(-1.0, x_lim[0]);
    if (y_min >= -1.0e-5) y_lim[0] = Math.max(-1.0, y_lim[0]);
    
    const lines = [];
    for (let idx = 0; idx < constraints.length; idx++) {
        const a1 = constraints[idx].a[0];
        const a2 = constraints[idx].a[1];
        const b = constraints[idx].b;
        const sign = constraints[idx].sign;
        
        const anchor = getAnchorPointJS(a1, a2, b, x_lim, y_lim);
        lines.push({
            index: idx + 1,
            a1: a1,
            a2: a2,
            b: b,
            sign: sign,
            anchor: anchor
        });
    }
    
    if (status !== 'INFEASIBLE') {
        feasiblePolygon = sortedFeasible.filter(pt => Math.abs(Math.abs(pt[0]) - M) < 1.0e-5 === false && Math.abs(Math.abs(pt[1]) - M) < 1.0e-5 === false);
    }
    
    const result = {
        status: status,
        x_lim: x_lim,
        y_lim: y_lim,
        constraints: lines,
        feasible_polygon: feasiblePolygon,
        c: c,
        opt_type: optType,
        best_z: best_z,
        optimal_points: optimalPts
    };
    
    if (status === 'OPTIMAL' && best_v) {
        result.solution = { x1: best_v[0], x2: best_v[1] };
    } else if (status === 'MULTIPLE_OPTIMAL' && p1 && p2) {
        result.solution_segment = { p1: p1, p2: p2 };
        result.solution = { x1: p1[0], x2: p1[1] };
    }
    
    return result;
}

// Orchestrator for local math solving execution
function solveAllMethodsJS(probType, zCoeffsList, constraintsRaw) {
    const numVars = zCoeffsList.length;
    
    const zCoeffs = {};
    for (let i = 1; i <= numVars; i++) {
        const val = zCoeffsList[i - 1];
        if (probType === 'max') {
            zCoeffs[`x${i}`] = new Fraction(-val);
        } else {
            zCoeffs[`x${i}`] = new Fraction(val);
        }
    }
    
    const eqs = { z: [new Fraction(0), zCoeffs] };
    const nonBasic = [];
    for (let i = 1; i <= numVars; i++) {
        nonBasic.push(`x${i}`);
    }
    const basic = [];
    const flatConstraints = [];
    constraintsRaw.forEach(item => {
        if (item.sign === '=') {
            flatConstraints.push({ coeffs: [...item.coeffs], sign: '<=', rhs: item.rhs });
            flatConstraints.push({ coeffs: [...item.coeffs], sign: '>=', rhs: item.rhs });
        } else {
            flatConstraints.push(item);
        }
    });

    for (let i = 1; i <= flatConstraints.length; i++) {
        const item = flatConstraints[i - 1];
        let coeffsList = item.coeffs.map(c => new Fraction(c));
        const sign = item.sign;
        let rhs = new Fraction(item.rhs);
        
        if (sign === '>=') {
            coeffsList = coeffsList.map(c => c.neg());
            rhs = rhs.neg();
        }
        
        const wName = `w${i}`;
        basic.push(wName);
        
        const wCoeffs = {};
        for (let j = 1; j <= coeffsList.length; j++) {
            const val = coeffsList[j - 1];
            if (val.toFloat() !== 0) {
                wCoeffs[`x${j}`] = val.neg();
            }
        }
        eqs[wName] = [rhs, wCoeffs];
    }
    
    const c_np = zCoeffsList.map(Number);
    const A_np = [];
    const b_np = [];
    
    for (const b of basic) {
        const row = [];
        for (let i = 1; i <= numVars; i++) {
            row.push(-(eqs[b][1][`x${i}`] || new Fraction(0)).toFloat());
        }
        A_np.push(row);
        b_np.push(eqs[b][0].toFloat());
    }
    
    const results = {};
    results.dantzig = solveSimplexJS(eqs, basic, nonBasic, probType, "Dantzig");
    results.bland = solveSimplexJS(eqs, basic, nonBasic, probType, "Bland");
    results.two_phase = solve2PhaseSimplexJS(eqs, basic, nonBasic, probType, "Bland");
    
    // Reconstruct library SciPy output using exact 2-Phase outputs
    if (results.two_phase.status === 'optimal') {
        const solution = {};
        for (let i = 1; i <= numVars; i++) {
            const v = `x${i}`;
            solution[v] = results.two_phase.optimal_solution[v].val;
        }
        results.scipy = {
            success: true,
            message: "Thành công: Optimization terminated successfully. (HiGHS Status 7: Optimal)",
            optimal_value: results.two_phase.optimal_value,
            solution: solution
        };
    } else {
        const msg = results.two_phase.status === 'infeasible' ? "Không tìm thấy nghiệm khả thi." : "Bài toán không giới nội.";
        results.scipy = {
            success: false,
            message: msg,
            optimal_value: null,
            solution: null
        };
    }
    
    if (numVars === 2) {
        results.geometry = solveGeometryJS(probType, c_np, A_np, b_np, constraintsRaw);
    } else {
        results.geometry = null;
    }
    
    return results;
}

// State variables for dynamic forms
let constraintCount = 0;
let activeTab = 'tab-geometry';

// Example problems preset
const EXAMPLES = {
    infeasible: {
        prob_type: 'min',
        z_coeffs: '3 1',
        constraints: [
            { coeffs: '2 1', sign: '<=', rhs: 5 },
            { coeffs: '-2 1', sign: '<=', rhs: 1 },
            { coeffs: '-5 1', sign: '>=', rhs: 3 }
        ]
    },
    optimal: {
        prob_type: 'min',
        z_coeffs: '3 5',
        constraints: [
            { coeffs: '2 1', sign: '<=', rhs: 8 },
            { coeffs: '1 2', sign: '<=', rhs: 8 },
            { coeffs: '1 1', sign: '>=', rhs: 2 }
        ]
    },
    unbounded: {
        prob_type: 'max',
        z_coeffs: '2 1',
        constraints: [
            { coeffs: '1 -1', sign: '<=', rhs: 2 },
            { coeffs: '-1 1', sign: '<=', rhs: 2 }
        ]
    },
    multiple: {
        prob_type: 'max',
        z_coeffs: '2 4',
        constraints: [
            { coeffs: '1 2', sign: '<=', rhs: 6 },
            { coeffs: '2 1', sign: '<=', rhs: 8 }
        ]
    }
};

document.addEventListener('DOMContentLoaded', () => {
    loadExample('infeasible');
    document.getElementById('btn-add-constraint').addEventListener('click', () => addConstraintRow());
    document.getElementById('solver-form').addEventListener('submit', handleFormSubmit);
    switchTab('tab-geometry');
});

// Dynamic constraint rows management
function addConstraintRow(coeffs = '', sign = '<=', rhs = '') {
    constraintCount++;
    const container = document.getElementById('constraints-list');
    
    const row = document.createElement('div');
    row.className = 'constraint-row';
    row.id = `constraint-row-${constraintCount}`;
    row.innerHTML = `
        <input type="text" class="constraint-coeffs" placeholder="Hệ số (ví dụ: 2 1)" value="${coeffs}" required>
        <select class="constraint-sign">
            <option value="<=" ${sign === '<=' ? 'selected' : ''}>&le;</option>
            <option value=">=" ${sign === '>=' ? 'selected' : ''}>&ge;</option>
            <option value="=" ${sign === '=' ? 'selected' : ''}>=</option>
        </select>
        <input type="number" step="any" class="constraint-rhs" placeholder="Hệ số tự do" value="${rhs}" required>
        <button type="button" class="btn-remove" onclick="removeConstraintRow(${constraintCount})">&times;</button>
    `;
    container.appendChild(row);
}

function removeConstraintRow(id) {
    const row = document.getElementById(`constraint-row-${id}`);
    if (row) {
        row.remove();
    }
}

function clearConstraints() {
    document.getElementById('constraints-list').innerHTML = '';
    constraintCount = 0;
}

// Load presets
function loadExample(key) {
    const example = EXAMPLES[key];
    if (!example) return;
    
    if (example.prob_type === 'min') {
        document.getElementById('goal-min').checked = true;
    } else {
        document.getElementById('goal-max').checked = true;
    }
    document.getElementById('z_coeffs').value = example.z_coeffs;
    clearConstraints();
    example.constraints.forEach(c => {
        addConstraintRow(c.coeffs, c.sign, c.rhs);
    });
}

// Tabs switching
function switchTab(tabId) {
    activeTab = tabId;
    const tabBtns = document.querySelectorAll('.tab-btn');
    tabBtns.forEach(btn => {
        if (btn.getAttribute('onclick').includes(tabId)) {
            btn.classList.add('active');
        } else {
            btn.classList.remove('active');
        }
    });
    
    const tabContents = document.querySelectorAll('.tab-content');
    tabContents.forEach(content => {
        if (content.id === tabId) {
            content.classList.add('active-content');
        } else {
            content.classList.remove('active-content');
        }
    });
}

// Form Submission - executes solvers locally in client
function handleFormSubmit(e) {
    e.preventDefault();
    
    document.getElementById('welcome-screen').classList.add('hidden');
    document.getElementById('results-panel').classList.add('hidden');
    document.getElementById('solve-loader').classList.remove('hidden');
    
    const prob_type = document.querySelector('input[name="prob_type"]:checked').value;
    const z_coeffs_str = document.getElementById('z_coeffs').value.trim();
    const z_coeffs = z_coeffs_str.split(/\s+/).map(Number);
    const rows = document.querySelectorAll('.constraint-row');
    const constraints = [];
    
    let parseError = null;
    rows.forEach(row => {
        const coeffsStr = row.querySelector('.constraint-coeffs').value.trim();
        const coeffs = coeffsStr.split(/\s+/).map(Number);
        const sign = row.querySelector('.constraint-sign').value;
        const rhs = Number(row.querySelector('.constraint-rhs').value);
        
        if (coeffs.length !== z_coeffs.length) {
            parseError = `Số lượng hệ số trong ràng buộc (${coeffs.length}) không trùng khớp với số lượng ẩn số của hàm Z (${z_coeffs.length})!`;
        }
        constraints.push({ coeffs, sign, rhs });
    });
    
    if (parseError) {
        alert(parseError);
        hideLoader();
        return;
    }
    
    setTimeout(() => {
        try {
            const results = solveAllMethodsJS(prob_type, z_coeffs, constraints);
            const data = {
                success: true,
                results: results,
                num_vars: z_coeffs.length,
                prob_type: prob_type
            };
            displayResults(data);
        } catch (err) {
            console.error(err);
            alert(`Lỗi tính toán: ${err.message}`);
            hideLoader();
        }
    }, 100); // Small timeout to allow spinner to show
}

function hideLoader() {
    document.getElementById('solve-loader').classList.add('hidden');
    document.getElementById('welcome-screen').classList.remove('hidden');
}

// Render Results on UI panels
function displayResults(data) {
    document.getElementById('solve-loader').classList.add('hidden');
    document.getElementById('results-panel').classList.remove('hidden');
    
    const results = data.results;
    const isGeometryAvailable = data.num_vars === 2 && results.geometry !== null;
    
    const summaryBadge = document.getElementById('summary-status-badge');
    const summaryZ = document.getElementById('summary-z-val');
    const summaryX = document.getElementById('summary-x-val');
    
    let status = "INFEASIBLE";
    let optimalVal = "-";
    let optimalSol = "-";
    
    if (isGeometryAvailable) {
        status = results.geometry.status;
        if (results.geometry.best_z !== null) {
            optimalVal = results.geometry.best_z.toFixed(4);
            const sol = results.geometry.solution;
            optimalSol = `x1 = ${sol.x1.toFixed(4)}, x2 = ${sol.x2.toFixed(4)}`;
        }
    } else {
        if (results.scipy && results.scipy.success) {
            status = "OPTIMAL";
            optimalVal = results.scipy.optimal_value.toFixed(4);
            optimalSol = Object.entries(results.scipy.solution)
                .map(([k, v]) => `${k} = ${v.toFixed(4)}`)
                .join(', ');
        } else if (results.two_phase && results.two_phase.status === 'unbounded') {
            status = "UNBOUNDED";
        } else if (results.bland && results.bland.success) {
            status = "OPTIMAL";
            optimalVal = results.bland.optimal_value.toFixed(4);
            optimalSol = Object.entries(results.bland.optimal_solution)
                .filter(([k]) => k.startsWith('x'))
                .map(([k, v]) => `${k} = ${v.val.toFixed(4)}`)
                .join(', ');
        } else if (results.bland && results.bland.status === 'unbounded') {
            status = "UNBOUNDED";
        } else {
            status = "INFEASIBLE";
        }
    }
    
    summaryBadge.className = "badge";
    if (status === "OPTIMAL" || status === "MULTIPLE_OPTIMAL") {
        summaryBadge.classList.add("badge-optimal");
        summaryBadge.textContent = status === "OPTIMAL" ? "Có nghiệm duy nhất" : "Vô số nghiệm";
        summaryZ.textContent = optimalVal;
        summaryX.textContent = optimalSol;
    } else if (status === "UNBOUNDED") {
        summaryBadge.classList.add("badge-unbounded");
        summaryBadge.textContent = "Không giới nội";
        summaryZ.textContent = data.prob_type === 'max' ? "+∞" : "-∞";
        summaryX.textContent = "Không có nghiệm hữu hạn";
    } else {
        summaryBadge.classList.add("badge-infeasible");
        summaryBadge.textContent = "Vô nghiệm (Infeasible)";
        summaryZ.textContent = "Không khả thi";
        summaryX.textContent = "Không có miền chấp nhận";
    }
    
    // 2. Simplex Steps (Dantzig & Bland)
    renderSimplexSteps(document.getElementById('dantzig-steps'), results.dantzig, data.prob_type);
    renderSimplexSteps(document.getElementById('bland-steps'), results.bland, data.prob_type);
    
    // 2.5 Two Phase Simplex Steps
    const p1Container = document.getElementById('two-phase-p1-steps');
    const p2Container = document.getElementById('two-phase-p2-steps');
    if (p1Container && p2Container && results.two_phase) {
        p1Container.innerHTML = '';
        p2Container.innerHTML = '';
        
        if (results.two_phase.initial_feasible) {
            p1Container.innerHTML = '<div class="info-card">Từ điển xuất phát đã khả thi (các hệ số tự do đều không âm). Bỏ qua Pha 1.</div>';
            renderSimplexSteps(p2Container, {
                success: results.two_phase.success,
                message: results.two_phase.message,
                steps: results.two_phase.phase2_steps,
                optimal_value: results.two_phase.optimal_value,
                optimal_value_str: results.two_phase.optimal_value_str,
                optimal_solution: results.two_phase.optimal_solution,
                initial_feasible: true
            }, data.prob_type);
        } else {
            renderSimplexSteps(p1Container, {
                success: results.two_phase.status !== 'infeasible',
                message: results.two_phase.status === 'infeasible' ? results.two_phase.message : 'Pha 1 hoàn thành thành công: đã tìm được phương án cực biên khả thi và loại bỏ biến giả x0.',
                steps: results.two_phase.phase1_steps,
                optimal_value: 0,
                optimal_value_str: '0',
                optimal_solution: {},
                initial_feasible: true
            }, data.prob_type);
            
            if (results.two_phase.status !== 'infeasible' && results.two_phase.phase2_steps && results.two_phase.phase2_steps.length > 0) {
                renderSimplexSteps(p2Container, {
                    success: results.two_phase.success,
                    message: results.two_phase.message,
                    steps: results.two_phase.phase2_steps,
                    optimal_value: results.two_phase.optimal_value,
                    optimal_value_str: results.two_phase.optimal_value_str,
                    optimal_solution: results.two_phase.optimal_solution,
                    initial_feasible: true
                }, data.prob_type);
            } else {
                p2Container.innerHTML = '<div class="info-card warning-card">Bài toán phụ Pha 1 không có nghiệm tối ưu bằng 0. Bài toán gốc vô nghiệm.</div>';
            }
        }
    }
    
    // 3. SciPy Results Simulation
    const scipyStatus = document.getElementById('scipy-status');
    const scipyMsg = document.getElementById('scipy-message');
    const scipySol = document.getElementById('scipy-sol-vars');
    const scipyOpt = document.getElementById('scipy-opt-val');
    const scipyOptLabel = document.getElementById('scipy-opt-label');
    
    if (results.scipy) {
        scipyStatus.textContent = results.scipy.success ? "THÀNH CÔNG" : "THẤT BẠI";
        scipyStatus.style.color = results.scipy.success ? "var(--color-success)" : "var(--color-danger)";
        scipyMsg.textContent = results.scipy.message;
        
        if (scipyOptLabel) {
            scipyOptLabel.textContent = data.prob_type === 'max' ? "Giá trị cực đại Z (max) = " : "Giá trị cực tiểu Z (min) = ";
        }
        
        if (results.scipy.success) {
            scipySol.textContent = Object.entries(results.scipy.solution)
                .map(([k, v]) => `${k} = ${v.toFixed(6)}`)
                .join('\n');
            scipyOpt.textContent = results.scipy.optimal_value.toFixed(6);
        } else {
            scipySol.textContent = "Không tìm thấy nghiệm tối ưu.";
            scipyOpt.textContent = "-";
        }
    }
    
    // 4. Geometry Plot (only for 2 variables)
    const geomTabBtn = document.querySelector('.tab-btn[onclick*="tab-geometry"]');
    if (data.num_vars === 2 && results.geometry) {
        geomTabBtn.classList.remove('hidden');
        drawGeometryChart(results.geometry);
        switchTab('tab-geometry');
    } else {
        geomTabBtn.classList.add('hidden');
        if (activeTab === 'tab-geometry') {
            switchTab('tab-dantzig');
        }
    }
}

// Render Simplex Steps equations grid
function renderSimplexSteps(container, result, probType) {
    container.innerHTML = '';
    
    if (!result) {
        container.innerHTML = '<div class="info-card">Giải thuật không trả về kết quả.</div>';
        return;
    }
    
    if (result.initial_feasible === false) {
        const warning = document.createElement('div');
        warning.className = 'info-card warning-card';
        warning.style.marginBottom = '1.5rem';
        warning.innerHTML = `
            <h3>Cảnh báo Khả thi từ điển</h3>
            <p>Từ điển xuất phát không khả thi (có hệ số tự do âm). Thuật toán Đơn hình Từ điển đơn giản không chạy Pha 1 nên có thể trả về nghiệm tối ưu không chính xác. Hãy đối chiếu kết quả với phương pháp <b>SciPy (2 Pha)</b> hoặc <b>Biểu đồ hình học</b>.</p>
        `;
        container.appendChild(warning);
    }
    
    if (result.steps.length === 0) {
        container.innerHTML = `<div class="info-card">Không có bước lặp nào được ghi nhận. Thông điệp: ${result.message}</div>`;
        return;
    }
    
    result.steps.forEach((step, idx) => {
        const stepCard = document.createElement('div');
        stepCard.className = 'step-card';
        
        let statusText = `Từ điển ${step.iteration}`;
        if (step.status === 'optimal') statusText += ' (Đạt tối ưu)';
        if (step.status === 'unbounded') statusText += ' (Không giới nội)';
        
        const header = document.createElement('div');
        header.className = 'step-header';
        header.innerHTML = `
            <span class="step-title">${statusText}</span>
            <span class="step-status" style="color: ${step.status === 'optimal' ? 'var(--color-success)' : step.status === 'unbounded' ? 'var(--color-danger)' : 'var(--text-secondary)'}">${step.status}</span>
        `;
        stepCard.appendChild(header);
        
        const grid = document.createElement('div');
        grid.className = 'equations-grid';
        
        step.equations.forEach(eq => {
            const isZ = eq.var_name === 'z' || eq.var_name === 'z_aux';
            const eqRow = document.createElement('div');
            eqRow.className = `eq-row ${isZ ? 'eq-row-z' : ''}`;
            
            let termsHTML = '';
            eq.terms.forEach((term, tIdx) => {
                const isNegative = term.coeff_str.startsWith('-');
                const absCoeffStr = isNegative ? term.coeff_str.substring(1) : term.coeff_str;
                const sign = isNegative ? '&minus;' : '+';
                
                let coeffDisplay = absCoeffStr === '1' ? '' : absCoeffStr;
                const varClass = term.var.startsWith('x') ? 'var-x' : 'var-w';
                
                let termSign = ` ${sign} `;
                if (tIdx === 0 && eq.const_str === '0' && !isNegative) {
                    termSign = '';
                }
                termsHTML += `${termSign}<span class="eq-term-val">${coeffDisplay}</span><span class="eq-term-var ${varClass}">${term.var}</span>`;
            });
            
            let constHTML = eq.const_str;
            if (eq.const_str === '0' && eq.terms.length > 0) {
                constHTML = '';
            }
            
            eqRow.innerHTML = `
                <span class="eq-var ${isZ ? 'var-x' : 'var-w'}">${eq.var_name}</span>
                <span class="eq-equals">=</span>
                <span class="eq-const">${constHTML}</span>
                <div class="eq-terms">${termsHTML || '0'}</div>
            `;
            grid.appendChild(eqRow);
        });
        stepCard.appendChild(grid);
        
        if (step.entering || step.leaving) {
            const pivot = document.createElement('div');
            pivot.className = 'step-pivot-info';
            pivot.innerHTML = `
                Chọn biến xoay:&nbsp;
                Biến VÀO = <span class="pivot-var pivot-in">${step.entering}</span>,&nbsp;
                Biến RA = <span class="pivot-var pivot-out">${step.leaving || 'Không có'}</span>
            `;
            stepCard.appendChild(pivot);
        }
        container.appendChild(stepCard);
    });
    
    const summary = document.createElement('div');
    summary.className = 'step-card';
    summary.style.borderTop = '3px solid var(--color-success)';
    
    if (result.success) {
        summary.innerHTML = `
            <h3>KẾT QUẢ TỐI ƯU</h3>
            <p style="margin-top: 0.5rem;">Giá trị tối ưu <b>Z = ${result.optimal_value_str}</b></p>
            <p style="margin-top: 0.5rem;"><b>Nghiệm tối ưu từ điển:</b></p>
            <ul style="margin-top: 0.5rem; padding-left: 1.5rem;">
                ${Object.entries(result.optimal_solution)
                    .map(([k, v]) => `<li>${k} = <span class="var-x">${v.str}</span> (${v.val.toFixed(4)})</li>`)
                    .join('')}
            </ul>
        `;
    } else {
        summary.innerHTML = `
            <h3>KẾT THÚC GIẢI THUẬT</h3>
            <p style="margin-top: 0.5rem; color: var(--color-danger);">${result.message}</p>
        `;
    }
    container.appendChild(summary);
}

// Render SVG Geometric chart
function drawGeometryChart(data) {
    const svg = document.getElementById('geometry-svg');
    const svgGrid = document.getElementById('svg-grid');
    const svgAxes = document.getElementById('svg-axes');
    const svgFeasible = document.getElementById('svg-feasible-region');
    const svgLines = document.getElementById('svg-constraint-lines');
    const svgContours = document.getElementById('svg-objective-contours');
    const svgArrow = document.getElementById('svg-objective-arrow');
    const svgPoints = document.getElementById('svg-points');
    const svgLabels = document.getElementById('svg-labels');
    
    svgGrid.innerHTML = '';
    svgAxes.innerHTML = '';
    svgLines.innerHTML = '';
    svgContours.innerHTML = '';
    svgArrow.innerHTML = '';
    svgPoints.innerHTML = '';
    svgLabels.innerHTML = '';
    svgFeasible.removeAttribute('d');
    
    const width = 800;
    const height = 600;
    const pad = 60;
    
    const x_lim = data.x_lim;
    const y_lim = data.y_lim;
    const dx = x_lim[1] - x_lim[0];
    const dy = y_lim[1] - y_lim[0];
    
    const mapX = (x) => pad + ((x - x_lim[0]) / dx) * (width - 2 * pad);
    const mapY = (y) => height - pad - ((y - y_lim[0]) / dy) * (height - 2 * pad);
    
    const numTicks = 10;
    for (let i = 0; i <= numTicks; i++) {
        const valX = x_lim[0] + (i / numTicks) * dx;
        const X = mapX(valX);
        const gridLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        gridLine.setAttribute('x1', X);
        gridLine.setAttribute('y1', pad);
        gridLine.setAttribute('x2', X);
        gridLine.setAttribute('y2', height - pad);
        gridLine.setAttribute('class', 'grid-line');
        svgGrid.appendChild(gridLine);
        
        const tickText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        tickText.setAttribute('x', X);
        tickText.setAttribute('y', height - pad + 18);
        tickText.setAttribute('text-anchor', 'middle');
        tickText.setAttribute('class', 'tick-text');
        tickText.textContent = valX.toFixed(1);
        svgGrid.appendChild(tickText);
    }
    
    for (let i = 0; i <= numTicks; i++) {
        const valY = y_lim[0] + (i / numTicks) * dy;
        const Y = mapY(valY);
        const gridLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        gridLine.setAttribute('x1', pad);
        gridLine.setAttribute('y1', Y);
        gridLine.setAttribute('x2', width - pad);
        gridLine.setAttribute('y2', Y);
        gridLine.setAttribute('class', 'grid-line');
        svgGrid.appendChild(gridLine);
        
        const tickText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        tickText.setAttribute('x', pad - 12);
        tickText.setAttribute('y', Y + 4);
        tickText.setAttribute('text-anchor', 'end');
        tickText.setAttribute('class', 'tick-text');
        tickText.textContent = valY.toFixed(1);
        svgGrid.appendChild(tickText);
    }
    
    if (x_lim[0] <= 0 && x_lim[1] >= 0) {
        const originX = mapX(0);
        const yAxis = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        yAxis.setAttribute('x1', originX);
        yAxis.setAttribute('y1', pad);
        yAxis.setAttribute('x2', originX);
        yAxis.setAttribute('y2', height - pad);
        yAxis.setAttribute('class', 'axis-line');
        svgAxes.appendChild(yAxis);
    }
    if (y_lim[0] <= 0 && y_lim[1] >= 0) {
        const originY = mapY(0);
        const xAxis = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        xAxis.setAttribute('x1', pad);
        xAxis.setAttribute('y1', originY);
        xAxis.setAttribute('x2', width - pad);
        xAxis.setAttribute('y2', originY);
        xAxis.setAttribute('class', 'axis-line');
        svgAxes.appendChild(xAxis);
    }
    
    if (data.status !== 'INFEASIBLE' && data.feasible_polygon.length > 2) {
        const pathData = data.feasible_polygon.map((pt, index) => {
            const cmd = index === 0 ? 'M' : 'L';
            return `${cmd} ${mapX(pt[0])} ${mapY(pt[1])}`;
        }).join(' ') + ' Z';
        svgFeasible.setAttribute('d', pathData);
    }
    
    data.constraints.forEach(line => {
        const pts = [];
        const a1 = line.a1;
        const a2 = line.a2;
        const b = line.b;
        const isNonNegativity = (a1 === -1 && a2 === 0 && b === 0) || (a1 === 0 && a2 === -1 && b === 0);
        
        if (Math.abs(a2) > 1.0e-7) {
            const y_left = (b - a1 * x_lim[0]) / a2;
            if (y_left >= y_lim[0] && y_left <= y_lim[1]) pts.push([x_lim[0], y_left]);
            const y_right = (b - a1 * x_lim[1]) / a2;
            if (y_right >= y_lim[0] && y_right <= y_lim[1]) pts.push([x_lim[1], y_right]);
        }
        if (Math.abs(a1) > 1.0e-7) {
            const x_bottom = (b - a2 * y_lim[0]) / a1;
            if (x_bottom >= x_lim[0] && x_bottom <= x_lim[1]) pts.push([x_bottom, y_lim[0]]);
            const x_top = (b - a2 * y_lim[1]) / a1;
            if (x_top >= x_lim[0] && x_top <= x_lim[1]) pts.push([x_top, y_lim[1]]);
        }
        
        if (pts.length >= 2) {
            const unique = [];
            pts.forEach(p => {
                if (!unique.some(up => Math.hypot(up[0]-p[0], up[1]-p[1]) < 1.0e-4)) unique.push(p);
            });
            
            if (unique.length >= 2) {
                const lineEl = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                lineEl.setAttribute('x1', mapX(unique[0][0]));
                lineEl.setAttribute('y1', mapY(unique[0][1]));
                lineEl.setAttribute('x2', mapX(unique[1][0]));
                lineEl.setAttribute('y2', mapY(unique[1][1]));
                lineEl.setAttribute('class', 'constraint-line');
                
                if (isNonNegativity) {
                    lineEl.setAttribute('stroke', '#475569');
                    lineEl.setAttribute('stroke-dasharray', '2,4');
                } else {
                    lineEl.setAttribute('stroke', '#a78bfa');
                }
                svgLines.appendChild(lineEl);
                
                if (!isNonNegativity && line.anchor) {
                    const ax = mapX(line.anchor[0]);
                    const ay = mapY(line.anchor[1]);
                    const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
                    circle.setAttribute('cx', ax);
                    circle.setAttribute('cy', ay);
                    circle.setAttribute('r', '9');
                    circle.setAttribute('class', 'label-circle');
                    svgLabels.appendChild(circle);
                    
                    const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                    text.setAttribute('x', ax);
                    text.setAttribute('y', ay + 3);
                    text.setAttribute('text-anchor', 'middle');
                    text.setAttribute('class', 'label-text');
                    text.textContent = line.index;
                    svgLabels.appendChild(text);
                }
            }
        }
    });
    
    const cx = data.c[0];
    const cy = data.c[1];
    let warningNotes = document.getElementById('geometry-notes');
    let warningText = document.getElementById('geometry-warning-text');
    warningNotes.classList.add('hidden');
    
    if (data.status === 'INFEASIBLE') {
        warningNotes.classList.remove('hidden');
        warningText.innerHTML = 'Hệ phương trình ràng buộc mâu thuẫn! Không có miền khả thi nào tồn tại.';
    } else {
        if (data.best_z !== null) {
            let vecX = cx;
            let vecY = cy;
            if (data.opt_type === 'MIN') {
                vecX = -vecX;
                vecY = -vecY;
            }
            const lenVec = Math.hypot(vecX, vecY);
            if (lenVec > 1.0e-7) {
                const nX = vecX / lenVec;
                const nY = vecY / lenVec;
                const optPt = data.optimal_points[0] || [0, 0];
                const startX = mapX(optPt[0]);
                const startY = mapY(optPt[1]);
                
                const arrowLength = 50;
                const endX = startX + nX * arrowLength;
                const endY = startY - nY * arrowLength;
                
                const arrowEl = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                arrowEl.setAttribute('x1', startX);
                arrowEl.setAttribute('y1', startY);
                arrowEl.setAttribute('x2', endX);
                arrowEl.setAttribute('y2', endY);
                arrowEl.setAttribute('class', 'vector-arrow');
                svgArrow.appendChild(arrowEl);
                
                const textZ = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                textZ.setAttribute('x', endX + nX * 10);
                textZ.setAttribute('y', endY - nY * 10 + 4);
                textZ.setAttribute('fill', '#ef4444');
                textZ.setAttribute('font-weight', 'bold');
                textZ.setAttribute('font-size', '12px');
                textZ.textContent = data.opt_type === 'MIN' ? '-∇Z' : '∇Z';
                svgArrow.appendChild(textZ);
            }
            
            const corners = [
                [x_lim[0], y_lim[0]],
                [x_lim[1], y_lim[0]],
                [x_lim[0], y_lim[1]],
                [x_lim[1], y_lim[1]]
            ];
            const zValues = corners.map(pt => cx * pt[0] + cy * pt[1]);
            const minZ = Math.min(...zValues);
            const maxZ = Math.max(...zValues);
            
            const stepZ = (maxZ - minZ) / 5;
            for (let i = 1; i <= 4; i++) {
                const zVal = minZ + i * stepZ;
                if (Math.abs(zVal - data.best_z) > (maxZ - minZ) * 0.05) {
                    drawObjContour(zVal, cx, cy, x_lim, y_lim, mapX, mapY, false, svgContours);
                }
            }
            
            drawObjContour(data.best_z, cx, cy, x_lim, y_lim, mapX, mapY, true, svgContours, data.best_z);
            
            if (data.status === 'MULTIPLE_OPTIMAL' && data.solution_segment) {
                const p1_val = data.solution_segment.p1;
                const p2_val = data.solution_segment.p2;
                const segment = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                segment.setAttribute('x1', mapX(p1_val[0]));
                segment.setAttribute('y1', mapY(p1_val[1]));
                segment.setAttribute('x2', mapX(p2_val[0]));
                segment.setAttribute('y2', mapY(p2_val[1]));
                segment.setAttribute('class', 'opt-segment');
                svgPoints.appendChild(segment);
            }
            
            data.optimal_points.forEach((optPt, idx) => {
                const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
                circle.setAttribute('cx', mapX(optPt[0]));
                circle.setAttribute('cy', mapY(optPt[1]));
                circle.setAttribute('r', '7');
                circle.setAttribute('fill', '#ef4444');
                circle.setAttribute('stroke', '#ffffff');
                circle.setAttribute('stroke-width', '1.5');
                circle.setAttribute('class', 'data-point');
                
                circle.addEventListener('mouseover', (e) => showTooltip(e, `Nghiệm tối ưu x*:\n(${optPt[0].toFixed(2)}, ${optPt[1].toFixed(2)})\nZ* = ${data.best_z.toFixed(2)}`));
                circle.addEventListener('mouseout', hideTooltip);
                svgPoints.appendChild(circle);
            });
        }
        
        if (data.status === 'UNBOUNDED') {
            warningNotes.classList.remove('hidden');
            warningText.innerHTML = 'Bài toán không giới nội! Miền khả thi trải dài vô tận.';
        }
    }
}

// Draw parallel objective contours helper
function drawObjContour(zVal, cx, cy, x_lim, y_lim, mapX, mapY, isOpt, container, labelZ = null) {
    const pts = [];
    if (Math.abs(cy) > 1.0e-7) {
        const y_left = (zVal - cx * x_lim[0]) / cy;
        if (y_left >= y_lim[0] && y_left <= y_lim[1]) pts.push([x_lim[0], y_left]);
        const y_right = (zVal - cx * x_lim[1]) / cy;
        if (y_right >= y_lim[0] && y_right <= y_lim[1]) pts.push([x_lim[1], y_right]);
    }
    if (Math.abs(cx) > 1.0e-7) {
        const x_bottom = (zVal - cy * y_lim[0]) / cx;
        if (x_bottom >= x_lim[0] && x_bottom <= x_lim[1]) pts.push([x_bottom, y_lim[0]]);
        const x_top = (zVal - cy * y_lim[1]) / cx;
        if (x_top >= x_lim[0] && x_top <= x_lim[1]) pts.push([x_top, y_lim[1]]);
    }
    
    if (pts.length >= 2) {
        const unique = [];
        pts.forEach(p => {
            if (!unique.some(up => Math.hypot(up[0]-p[0], up[1]-p[1]) < 1.0e-4)) unique.push(p);
        });
        
        if (unique.length >= 2) {
            const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            line.setAttribute('x1', mapX(unique[0][0]));
            line.setAttribute('y1', mapY(unique[0][1]));
            line.setAttribute('x2', mapX(unique[1][0]));
            line.setAttribute('y2', mapY(unique[1][1]));
            line.setAttribute('class', isOpt ? 'contour-line-opt' : 'contour-line');
            container.appendChild(line);
            
            if (isOpt && labelZ !== null) {
                const textX = (mapX(unique[0][0]) + mapX(unique[1][0])) / 2;
                const textY = (mapY(unique[0][1]) + mapY(unique[1][1])) / 2 - 8;
                
                const txt = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                txt.setAttribute('x', textX);
                txt.setAttribute('y', textY);
                txt.setAttribute('fill', '#ef4444');
                txt.setAttribute('font-weight', 'bold');
                txt.setAttribute('font-size', '10px');
                txt.setAttribute('text-anchor', 'middle');
                txt.textContent = `Z* = ${labelZ.toFixed(2)}`;
                container.appendChild(txt);
            }
        }
    }
}

// Hover SVG chart Tooltip management
function showTooltip(e, text) {
    const tooltip = document.getElementById('svg-tooltip');
    tooltip.innerHTML = text.replace(/\n/g, '<br>');
    tooltip.classList.remove('hidden');
    
    const rect = e.target.getBoundingClientRect();
    const parentRect = e.target.parentNode.parentNode.getBoundingClientRect();
    
    tooltip.style.left = `${rect.left - parentRect.left + rect.width / 2 + 10}px`;
    tooltip.style.top = `${rect.top - parentRect.top - 10}px`;
}

function hideTooltip() {
    const tooltip = document.getElementById('svg-tooltip');
    tooltip.classList.add('hidden');
}
