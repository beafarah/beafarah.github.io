---
layout: default
permalink: /power/
title: power calculator
nav: true
nav_order: 7
---

<h2>power for test of equality of quantiles</h2>

Power calculator based on the test of equality of quantiles as described in the paper ["Univariate and multivariate tests of equality of quantiles with right-censored data"](https://arxiv.org/abs/2505.03234).
This calculator assumes exponential distributions for control and censoring times, with balanced groups. The default parameters are the ones that correspond to some of the simulations from section 3.1 of our paper.

You can choose the model type: *Proportional hazards* (exponential experimental arm), and *Nonproportional hazards with late treatment effects* (piecewise exponential experimental arm).

Enter the parameters and click on **Calculate Power** to view the analytical power alongside the survival curves :)

<p><b>Choose calculation type:</b></p>
<button onclick="setCalcType('power')">Power Calculation</button>
<button onclick="setCalcType('samplesize')">Sample Size Calculation</button>
<p>Current calculation: <span id="calc-type">Power</span></p>

<p><b>Choose model type:</b></p>
<button onclick="setModel('exponential')">Proportional</button>
<button onclick="setModel('piecewise')">Nonproportional</button>
<p>Current model: <span id="model-type">Exponential</span></p>

<form id="power-form">
  <label>Probability: <input type="number" id="prob" step="any" required value="0.5"></label><br>
  <label>Difference: <input type="number" id="diff" step="any" required value="0.1"></label><br>

  <div id="sample-size-block">
    <label>Total Sample Size: <input type="number" id="sample-size" required value="1000"></label><br>
  </div>

  <div id="desired-power-block" style="display:none">
    <label>Desired Power: <input type="number" id="desired-power" step="any" required value="0.9"></label><br>
  </div>

  <label>Rate of Control Arm: <input type="number" id="rate-control" step="any" required value="1.5"></label><br>
  <label>Rate of Censoring: <input type="number" id="rate-cens" step="any" required value="0.48"></label><br>
  <label>Significance Level (alpha): <input type="number" id="alpha" step="any" required value="0.05"></label><br>

  <div id="piecewise-options" style="display:none">
    <label>Time Cutoff (tcut): <input type="number" id="tcut" step="any" value="0.2"></label><br>
  </div>

  <button type="submit">Calculate</button>
</form>

<p id="result"></p>

<canvas id="survival-chart" width="800" height="400"></canvas>

<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-annotation@1.4.0/dist/chartjs-plugin-annotation.min.js"></script>

{% raw %}
<script>

let model = 'exponential';
let calcType = 'power';

function setModel(m) {
  model = m;
  document.getElementById("model-type").innerText = m.charAt(0).toUpperCase() + m.slice(1);
  document.getElementById("piecewise-options").style.display = (m === 'piecewise') ? 'block' : 'none';
}

function setCalcType(t) {
  calcType = t;
  document.getElementById("calc-type").innerText = (t === 'power' ? "Power" : "Sample Size");
  document.getElementById("sample-size-block").style.display = (t === 'power') ? 'block' : 'none';
  document.getElementById("desired-power-block").style.display = (t === 'samplesize') ? 'block' : 'none';
}


/* ---- Utility functions ---- */

function normCDF(x) {
  var sign = x < 0 ? -1 : 1;
  x = Math.abs(x) / Math.sqrt(2);
  var a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741,
      a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
  var t = 1 / (1 + p * x);
  var y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
  return 0.5 * (1 + sign * y);
}

function inverseErf(x) {
  let a = 0.147;
  let ln = Math.log(1 - x * x);
  let term1 = 2 / (Math.PI * a) + ln / 2;
  let term2 = ln / a;
  return Math.sign(x) * Math.sqrt(Math.sqrt(term1 * term1 - term2) - term1);
}

function normSInv(p) {
  return Math.sqrt(2) * inverseErf(2 * p - 1);
}

function expo_pdf(x, lambda) {
  return lambda * Math.exp(-lambda * x);
}

function piecewise_pdf(x, rateC, rateE, tcut) {
  if (x <= tcut) return rateC * Math.exp(-rateC * x);
  const s_t = Math.exp(-rateC * tcut);
  return rateE * s_t * Math.exp(-rateE * (x - tcut));
}

/* ---- CORE POWER COMPUTATION (returns power and needed quantities) ---- */

function computePower(p, diff, rateC, rateCens, alpha, tcut, n) {

  const z_crit = Math.abs(normSInv(1 - alpha/2));
  const quantC = -Math.log(1 - p) / rateC;

  let rateE, quantE, phiE;

  if (model === 'exponential') {
    rateE = -Math.log(1 - p) / (quantC - diff);
    quantE = quantC - diff;
    phiE = rateE / (rateE + rateCens) * (Math.exp((rateE + rateCens)*quantE) - 1);
  } else {
    if (quantC - tcut <= diff) return null;
    const boundary = 1 - Math.exp(-rateC * tcut);

    rateE = (Math.log(1 - p) + rateC*tcut) / (tcut + diff - quantC);

    quantE = (p < boundary)
       ? -Math.log(1 - p) / rateC
       : tcut - ((Math.log(1 - p) + rateC*tcut) / rateE);

    phiE =
      (rateC/(rateC+rateCens))*(Math.exp((rateC+rateCens)*tcut)-1) +
      (rateE/(rateE+rateCens))*Math.exp((rateC-rateE)*tcut)*
      (Math.exp((rateE+rateCens)*quantE) - Math.exp((rateE+rateCens)*tcut));
  }

  const phiC = rateC/(rateC+rateCens) * (Math.exp((rateC+rateCens)*quantC)-1);

  let sigma2;
  if (model === 'exponential') {
    sigma2 = Math.pow(1 - p,2) *
       (phiC/(0.5*Math.pow(expo_pdf(quantC,rateC),2)) +
        phiE/(0.5*Math.pow(expo_pdf(quantE,rateE),2)));
  } else {
    sigma2 = Math.pow(1 - p,2) *
       (phiC/(0.5*Math.pow(expo_pdf(quantC,rateC),2)) +
        phiE/(0.5*Math.pow(piecewise_pdf(quantE,rateC,rateE,tcut),2)));
  }

  const se = Math.sqrt(sigma2/n);
  const power = 1 - normCDF(z_crit - diff/se) + normCDF(-z_crit - diff/se);

  return {power, quantC, quantE, rateE};
}


/* ---- MAIN FORM HANDLER ---- */

window.addEventListener("DOMContentLoaded", () => {
  const form = document.getElementById("power-form");

  form.addEventListener("submit", function(e) {
    e.preventDefault();

    const p = parseFloat(prob.value);
    const diff = parseFloat(document.getElementById("diff").value);
    const rateC = parseFloat(document.getElementById("rate-control").value);
    const rateCens = parseFloat(document.getElementById("rate-cens").value);
    const alpha = parseFloat(document.getElementById("alpha").value);
    const tcut = (model === 'piecewise') ? parseFloat(document.getElementById("tcut").value) : null;

    let result;

    /* ------------------- POWER CALCULATION MODE ------------------- */
    if (calcType === "power") {
      const n = parseFloat(document.getElementById("sample-size").value);
      result = computePower(p, diff, rateC, rateCens, alpha, tcut, n);
      if (!result) { alert("Invalid parameters"); return; }
      document.getElementById("result").innerText =
        `Estimated Power: ${(result.power*100).toFixed(2)}%`;
    }

    /* ---------------- SAMPLE SIZE CALCULATION MODE ---------------- */
    else {
      const desiredPower = parseFloat(document.getElementById("desired-power").value);

      let n = 10;  
      let pow = 0;

      while (pow < desiredPower && n < 500000) {
        let R = computePower(p, diff, rateC, rateCens, alpha, tcut, n);
        if (!R) { alert("Invalid parameters."); return; }
        pow = R.power;
        n += 2;   // increase resolution as you wish
      }

      document.getElementById("result").innerText =
        `Minimum sample size per arm = ${Math.ceil(n/2)} (total n = ${n}, achieved power ${(pow*100).toFixed(2)}%)`;

      result = computePower(p, diff, rateC, rateCens, alpha, tcut, n);
    }

    /* ---- Plot survival curves ---- */

    const quantC = result.quantC;
    const rateE = result.rateE;

    const timeMax = quantC * 1.5;
    const timePoints = Array.from({length:100}, (_,i)=> +(timeMax*i/99).toFixed(2));

    const survivalC = timePoints.map(t => Math.exp(-rateC*t));
    const survivalE = (model==="exponential")
         ? timePoints.map(t => Math.exp(-rateE*t))
         : timePoints.map(t => t<=tcut
              ? Math.exp(-rateC*t)
              : Math.exp(-rateC*tcut)*Math.exp(-rateE*(t-tcut)));

    const ctx = document.getElementById("survival-chart").getContext("2d");
    if (window.survivalChartInstance) window.survivalChartInstance.destroy();
    window.survivalChartInstance = new Chart(ctx, {
      type: "line",
      data: {
        labels: timePoints,
        datasets: [
          { label: "Control Arm", data: survivalC, borderColor: "limegreen", fill: false, tension: 0.3, borderWidth: 2 },
          { label: "Experimental Arm", data: survivalE, borderColor: "darkgreen", fill: false, tension: 0.3, borderWidth: 2 }
        ]
      },
      options: {
        responsive: true,
        plugins: {
          title: { display: true, text: "Survival Functions", font: { size: 18 } },
          legend: { labels: { font: { size: 14 } } },
          tooltip: {
            callbacks: {
              label: function(context) {
                      let yVal = context.raw; // use raw value
                return `${context.dataset.label}: ${yVal.toFixed(2)}`;
              }
            }
          },
          annotation: {
            annotations: {
              hLine: {
                type: 'line',
                yMin: 1 - p,
                yMax: 1 - p,
                borderColor: 'green',
                borderWidth: 2,
                borderDash: [6, 6],
                label: {
                  content: `1 - p = ${(1 - p).toFixed(2)}`,
                  enabled: true,
                  position: 'start',
                  backgroundColor: 'rgba(0,0,0,0.7)',
                  color: '#fff',
                  font: { style: 'italic' }
                }
              }
            }
          }
        },
        scales: {
          x: { title: { display: true, text: "Time", font: { size: 16 } },
             ticks: { stepSize: 0.2, callback: val => val.toFixed(1) }
             },
          y: {
            min: 0, max: 1,
            title: { display: true, text: "Survival Probability", font: { size: 16 } },
            ticks: { stepSize: 0.2, callback: val => val.toFixed(1) }
          }
        }
      }
    });
  });
});
</script>
{% endraw %}

