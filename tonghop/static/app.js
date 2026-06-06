// State variables
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
    // Add default constraints
    loadExample('infeasible'); // Start with the user request case
    
    // Bind buttons
    document.getElementById('btn-add-constraint').addEventListener('click', () => addConstraintRow());
    document.getElementById('solver-form').addEventListener('submit', handleFormSubmit);
    
    // Setup tabs
    switchTab('tab-geometry');
});

// Dynamic constraint rows
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
    
    // Set type
    if (example.prob_type === 'min') {
        document.getElementById('goal-min').checked = true;
    } else {
        document.getElementById('goal-max').checked = true;
    }
    
    // Set Z
    document.getElementById('z_coeffs').value = example.z_coeffs;
    
    // Set constraints
    clearConstraints();
    example.constraints.forEach(c => {
        addConstraintRow(c.coeffs, c.sign, c.rhs);
    });
}

// Tabs switching
function switchTab(tabId) {
    activeTab = tabId;
    
    // Update active tab button style
    const tabBtns = document.querySelectorAll('.tab-btn');
    tabBtns.forEach(btn => {
        if (btn.getAttribute('onclick').includes(tabId)) {
            btn.classList.add('active');
        } else {
            btn.classList.remove('active');
        }
    });
    
    // Update active content visibility
    const tabContents = document.querySelectorAll('.tab-content');
    tabContents.forEach(content => {
        if (content.id === tabId) {
            content.classList.add('active-content');
        } else {
            content.classList.remove('active-content');
        }
    });
}

// Form Submission & API call
async function handleFormSubmit(e) {
    e.preventDefault();
    
    // Show loader
    document.getElementById('welcome-screen').classList.add('hidden');
    document.getElementById('results-panel').classList.add('hidden');
    document.getElementById('solve-loader').classList.remove('hidden');
    
    // Collect data
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
    
    const payload = {
        prob_type,
        z_coeffs,
        constraints
    };
    
    try {
        const response = await fetch('/api/solve', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(payload)
        });
        
        const data = await response.json();
        if (data.success) {
            displayResults(data);
        } else {
            alert(`Lỗi: ${data.error}`);
            hideLoader();
        }
    } catch (err) {
        alert(`Không thể kết nối đến máy chủ: ${err}`);
        hideLoader();
    }
}

function hideLoader() {
    document.getElementById('solve-loader').classList.add('hidden');
    document.getElementById('welcome-screen').classList.remove('hidden');
}

// Display results on UI
function displayResults(data) {
    document.getElementById('solve-loader').classList.add('hidden');
    document.getElementById('results-panel').classList.remove('hidden');
    
    const results = data.results;
    const isGeometryAvailable = data.num_vars === 2 && results.geometry !== null;
    
    // 1. General Summary
    const summaryBadge = document.getElementById('summary-status-badge');
    const summaryZ = document.getElementById('summary-z-val');
    const summaryX = document.getElementById('summary-x-val');
    
    // Find representative status
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
        // Fallback to Scipy or Simplex
        if (results.scipy && results.scipy.success) {
            status = "OPTIMAL";
            optimalVal = results.scipy.optimal_value.toFixed(4);
            optimalSol = Object.entries(results.scipy.solution)
                .map(([k, v]) => `${k} = ${v.toFixed(4)}`)
                .join(', ');
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
    
    // Styling status badge
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
    
    // 3. SciPy Results
    const scipyStatus = document.getElementById('scipy-status');
    const scipyMsg = document.getElementById('scipy-message');
    const scipySol = document.getElementById('scipy-sol-vars');
    const scipyOpt = document.getElementById('scipy-opt-val');
    
    if (results.scipy) {
        scipyStatus.textContent = results.scipy.success ? "THÀNH CÔNG" : "THẤT BẠI";
        scipyStatus.style.color = results.scipy.success ? "var(--color-success)" : "var(--color-danger)";
        scipyMsg.textContent = results.scipy.message;
        
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

// Render Simplex steps into HTML table
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
        
        // Equations rendering
        const grid = document.createElement('div');
        grid.className = 'equations-grid';
        
        step.equations.forEach(eq => {
            const isZ = eq.var_name === 'z';
            const eqRow = document.createElement('div');
            eqRow.className = `eq-row ${isZ ? 'eq-row-z' : ''}`;
            
            // Format equation: var_name = const + term1 + term2...
            let termsHTML = '';
            eq.terms.forEach((term, tIdx) => {
                const isNegative = term.coeff_str.startsWith('-');
                const absCoeffStr = isNegative ? term.coeff_str.substring(1) : term.coeff_str;
                const sign = isNegative ? '&minus;' : '+';
                
                let coeffDisplay = absCoeffStr === '1' ? '' : absCoeffStr;
                const varClass = term.var.startsWith('x') ? 'var-x' : 'var-w';
                
                // For first term, if positive, skip sign representation if constant is 0
                let termSign = ` ${sign} `;
                if (tIdx === 0 && eq.const_str === '0' && !isNegative) {
                    termSign = '';
                }
                
                termsHTML += `${termSign}<span class="eq-term-val">${coeffDisplay}</span><span class="eq-term-var ${varClass}">${term.var}</span>`;
            });
            
            // Constant handling
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
        
        // Pivot variable information
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
    
    // Result summary at bottom of steps
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

// Draw premium SVG chart
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
    
    // Reset SVG elements
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
    const pad = 60; // Margins
    
    const x_lim = data.x_lim;
    const y_lim = data.y_lim;
    const dx = x_lim[1] - x_lim[0];
    const dy = y_lim[1] - y_lim[0];
    
    // Linear coordinates mapping
    const mapX = (x) => pad + ((x - x_lim[0]) / dx) * (width - 2 * pad);
    const mapY = (y) => height - pad - ((y - y_lim[0]) / dy) * (height - 2 * pad);
    
    // Invert mapping for hover coordinates
    const invX = (X) => x_lim[0] + ((X - pad) / (width - 2 * pad)) * dx;
    const invY = (Y) => y_lim[0] + ((height - pad - Y) / (height - 2 * pad)) * dy;
    
    // Draw grid lines
    const numTicks = 10;
    // X-grid & ticks
    for (let i = 0; i <= numTicks; i++) {
        const valX = x_lim[0] + (i / numTicks) * dx;
        const X = mapX(valX);
        
        // Vertical grid line
        const gridLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        gridLine.setAttribute('x1', X);
        gridLine.setAttribute('y1', pad);
        gridLine.setAttribute('x2', X);
        gridLine.setAttribute('y2', height - pad);
        gridLine.setAttribute('class', 'grid-line');
        svgGrid.appendChild(gridLine);
        
        // Label on axis
        const tickText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        tickText.setAttribute('x', X);
        tickText.setAttribute('y', height - pad + 18);
        tickText.setAttribute('text-anchor', 'middle');
        tickText.setAttribute('class', 'tick-text');
        tickText.textContent = valX.toFixed(1);
        svgGrid.appendChild(tickText);
    }
    
    // Y-grid & ticks
    for (let i = 0; i <= numTicks; i++) {
        const valY = y_lim[0] + (i / numTicks) * dy;
        const Y = mapY(valY);
        
        // Horizontal grid line
        const gridLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        gridLine.setAttribute('x1', pad);
        gridLine.setAttribute('y1', Y);
        gridLine.setAttribute('x2', width - pad);
        gridLine.setAttribute('y2', Y);
        gridLine.setAttribute('class', 'grid-line');
        svgGrid.appendChild(gridLine);
        
        // Label on axis
        const tickText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        tickText.setAttribute('x', pad - 12);
        tickText.setAttribute('y', Y + 4);
        tickText.setAttribute('text-anchor', 'end');
        tickText.setAttribute('class', 'tick-text');
        tickText.textContent = valY.toFixed(1);
        svgGrid.appendChild(tickText);
    }
    
    // Draw major axes (where x = 0 and y = 0 if they sit inside limits)
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
    
    // Draw feasible region polygon
    if (data.status !== 'INFEASIBLE' && data.feasible_polygon.length > 2) {
        const pathData = data.feasible_polygon.map((pt, index) => {
            const cmd = index === 0 ? 'M' : 'L';
            return `${cmd} ${mapX(pt[0])} ${mapY(pt[1])}`;
        }).join(' ') + ' Z';
        svgFeasible.setAttribute('d', pathData);
    }
    
    // Draw constraint lines
    data.constraints.forEach(line => {
        // Find line intersection with borders of limits
        const pts = [];
        const a1 = line.a1;
        const a2 = line.a2;
        const b = line.b;
        
        // Is standard constraint line or axis non-negativity boundary?
        // Note: we identify non-negativity constraints like x1 >= 0 (represented as -x1 <= 0)
        const isNonNegativity = (a1 === -1 && a2 === 0 && b === 0) || (a1 === 0 && a2 === -1 && b === 0);
        
        if (Math.abs(a2) > 1e-7) {
            // y = (b - a1*x)/a2
            const y_left = (b - a1 * x_lim[0]) / a2;
            if (y_left >= y_lim[0] && y_left <= y_lim[1]) pts.push([x_lim[0], y_left]);
            
            const y_right = (b - a1 * x_lim[1]) / a2;
            if (y_right >= y_lim[0] && y_right <= y_lim[1]) pts.push([x_lim[1], y_right]);
        }
        if (Math.abs(a1) > 1e-7) {
            // x = (b - a2*y)/a1
            const x_bottom = (b - a2 * y_lim[0]) / a1;
            if (x_bottom >= x_lim[0] && x_bottom <= x_lim[1]) pts.push([x_bottom, y_lim[0]]);
            
            const x_top = (b - a2 * y_lim[1]) / a1;
            if (x_top >= x_lim[0] && x_top <= x_lim[1]) pts.push([x_top, y_lim[1]]);
        }
        
        // Sort and draw the line segment inside viewbox
        if (pts.length >= 2) {
            // Unique points
            const unique = [];
            pts.forEach(p => {
                if (!unique.some(up => Math.hypot(up[0]-p[0], up[1]-p[1]) < 1e-4)) {
                    unique.push(p);
                }
            });
            
            if (unique.length >= 2) {
                const lineEl = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                lineEl.setAttribute('x1', mapX(unique[0][0]));
                lineEl.setAttribute('y1', mapY(unique[0][1]));
                lineEl.setAttribute('x2', mapX(unique[1][0]));
                lineEl.setAttribute('y2', mapY(unique[1][1]));
                lineEl.setAttribute('class', 'constraint-line');
                
                // Highlight different types of constraints
                if (isNonNegativity) {
                    lineEl.setAttribute('stroke', '#475569');
                    lineEl.setAttribute('stroke-dasharray', '2,4');
                } else {
                    lineEl.setAttribute('stroke', '#a78bfa'); // Violet accent for standard constraints
                }
                
                svgLines.appendChild(lineEl);
                
                // Label circle and text
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
    
    // Draw objective function contours and direction arrow
    const cx = data.c[0];
    const cy = data.c[1];
    
    let warningNotes = document.getElementById('geometry-notes');
    let warningText = document.getElementById('geometry-warning-text');
    warningNotes.classList.add('hidden');
    
    if (data.status === 'INFEASIBLE') {
        warningNotes.classList.remove('hidden');
        warningText.innerHTML = 'Hệ phương trình ràng buộc mâu thuẫn! Không có vùng giao giữa tất cả các ràng buộc. Không có nghiệm khả thi nào tồn tại.';
    } else {
        // Draw contour lines for Z
        if (data.best_z !== null) {
            // Objective direction vector: c vector
            let vecX = cx;
            let vecY = cy;
            if (data.opt_type === 'MIN') {
                vecX = -vecX;
                vecY = -vecY;
            }
            // Normalize direction vector
            const lenVec = Math.hypot(vecX, vecY);
            if (lenVec > 1e-7) {
                const nX = vecX / lenVec;
                const nY = vecY / lenVec;
                
                // Draw vector arrow starting from optimal point or origin
                const optPt = data.optimal_points[0] || [0, 0];
                const startX = mapX(optPt[0]);
                const startY = mapY(optPt[1]);
                
                // Arrow length (about 60px)
                const arrowLength = 50;
                const endX = startX + nX * arrowLength;
                const endY = startY - nY * arrowLength; // y-axis inverted
                
                const arrowEl = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                arrowEl.setAttribute('x1', startX);
                arrowEl.setAttribute('y1', startY);
                arrowEl.setAttribute('x2', endX);
                arrowEl.setAttribute('y2', endY);
                arrowEl.setAttribute('class', 'vector-arrow');
                svgArrow.appendChild(arrowEl);
                
                // Add "Z" text near arrow
                const textZ = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                textZ.setAttribute('x', endX + nX * 10);
                textZ.setAttribute('y', endY - nY * 10 + 4);
                textZ.setAttribute('fill', '#ef4444');
                textZ.setAttribute('font-weight', 'bold');
                textZ.setAttribute('font-size', '12px');
                textZ.textContent = data.opt_type === 'MIN' ? '-∇Z' : '∇Z';
                svgArrow.appendChild(textZ);
            }
            
            // Draw objective contours
            // We calculate Z at all four corners of viewbox to find Z min/max
            const corners = [
                [x_lim[0], y_lim[0]],
                [x_lim[1], y_lim[0]],
                [x_lim[0], y_lim[1]],
                [x_lim[1], y_lim[1]]
            ];
            const zValues = corners.map(pt => cx * pt[0] + cy * pt[1]);
            const minZ = Math.min(...zValues);
            const maxZ = Math.max(...zValues);
            
            // Plot 4 parallel objective lines
            const stepZ = (maxZ - minZ) / 5;
            for (let i = 1; i <= 4; i++) {
                const zVal = minZ + i * stepZ;
                // Avoid drawing extremely close to optimal Z to avoid clutter
                if (Math.abs(zVal - data.best_z) > (maxZ - minZ) * 0.05) {
                    drawObjContour(zVal, cx, cy, x_lim, y_lim, mapX, mapY, false, svgContours);
                }
            }
            
            // Draw OPTIMAL Z objective line
            drawObjContour(data.best_z, cx, cy, x_lim, y_lim, mapX, mapY, true, svgContours, data.best_z);
            
            // Highlight multiple optimal solutions segment
            if (data.status === 'MULTIPLE_OPTIMAL' && data.solution_segment) {
                const p1 = data.solution_segment.p1;
                const p2 = data.solution_segment.p2;
                const segment = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                segment.setAttribute('x1', mapX(p1[0]));
                segment.setAttribute('y1', mapY(p1[1]));
                segment.setAttribute('x2', mapX(p2[0]));
                segment.setAttribute('y2', mapY(p2[1]));
                segment.setAttribute('class', 'opt-segment');
                svgPoints.appendChild(segment);
            }
            
            // Draw optimal points as star
            data.optimal_points.forEach((optPt, idx) => {
                const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
                circle.setAttribute('cx', mapX(optPt[0]));
                circle.setAttribute('cy', mapY(optPt[1]));
                circle.setAttribute('r', '7');
                circle.setAttribute('fill', '#ef4444');
                circle.setAttribute('stroke', '#ffffff');
                circle.setAttribute('stroke-width', '1.5');
                circle.setAttribute('class', 'data-point');
                
                // Hover effect tooltips
                circle.addEventListener('mouseover', (e) => showTooltip(e, `Nghiệm tối ưu x*:\n(${optPt[0].toFixed(2)}, ${optPt[1].toFixed(2)})\nZ* = ${data.best_z.toFixed(2)}`));
                circle.addEventListener('mouseout', hideTooltip);
                
                svgPoints.appendChild(circle);
            });
        }
        
        if (data.status === 'UNBOUNDED') {
            warningNotes.classList.remove('hidden');
            warningText.innerHTML = 'Bài toán không giới nội! Miền chấp nhận trải dài vô tận theo hướng tăng/giảm hàm mục tiêu. Z tăng hoặc giảm về vô cùng.';
        }
    }
    
    // SVG tooltip binding
    svg.addEventListener('mousemove', (e) => {
        // Optional grid tracking details
    });
}

// Helper to draw objective contours
function drawObjContour(zVal, cx, cy, x_lim, y_lim, mapX, mapY, isOpt, container, labelZ = null) {
    const pts = [];
    if (Math.abs(cy) > 1e-7) {
        // y = (z - cx*x)/cy
        const y_left = (zVal - cx * x_lim[0]) / cy;
        if (y_left >= y_lim[0] && y_left <= y_lim[1]) pts.push([x_lim[0], y_left]);
        
        const y_right = (zVal - cx * x_lim[1]) / cy;
        if (y_right >= y_lim[0] && y_right <= y_lim[1]) pts.push([x_lim[1], y_right]);
    }
    if (Math.abs(cx) > 1e-7) {
        // x = (z - cy*y)/cx
        const x_bottom = (zVal - cy * y_lim[0]) / cx;
        if (x_bottom >= x_lim[0] && x_bottom <= x_lim[1]) pts.push([x_bottom, y_lim[0]]);
        
        const x_top = (zVal - cy * y_lim[1]) / cx;
        if (x_top >= x_lim[0] && x_top <= x_lim[1]) pts.push([x_top, y_lim[1]]);
    }
    
    if (pts.length >= 2) {
        const unique = [];
        pts.forEach(p => {
            if (!unique.some(up => Math.hypot(up[0]-p[0], up[1]-p[1]) < 1e-4)) {
                unique.push(p);
            }
        });
        
        if (unique.length >= 2) {
            const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            line.setAttribute('x1', mapX(unique[0][0]));
            line.setAttribute('y1', mapY(unique[0][1]));
            line.setAttribute('x2', mapX(unique[1][0]));
            line.setAttribute('y2', mapY(unique[1][1]));
            line.setAttribute('class', isOpt ? 'contour-line-opt' : 'contour-line');
            container.appendChild(line);
            
            // Text annotation for optimal contour
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

// Tooltip helpers
function showTooltip(e, text) {
    const tooltip = document.getElementById('svg-tooltip');
    tooltip.innerHTML = text.replace(/\n/g, '<br>');
    tooltip.classList.remove('hidden');
    
    // Position
    const rect = e.target.getBoundingClientRect();
    const parentRect = e.target.parentNode.parentNode.getBoundingClientRect();
    
    tooltip.style.left = `${rect.left - parentRect.left + rect.width / 2 + 10}px`;
    tooltip.style.top = `${rect.top - parentRect.top - 10}px`;
}

function hideTooltip() {
    const tooltip = document.getElementById('svg-tooltip');
    tooltip.classList.add('hidden');
}
