// ==UserScript==
// @name         SpectrumViwerLinkInjector
// @namespace    http://tampermonkey.net/
// @version      0.1
// @description  try to take over the world!
// @author       pgarrett@scripps.edu
// @match        http://ip2.scripps.edu/ip2/showLorikeet.html*
// @match        http://172.29.227.121/ip2/showLorikeet.html*
// @icon         https://www.google.com/s2/favicons?sz=64&domain=scripps.edu
// @require      https://cdnjs.cloudflare.com/ajax/libs/lz-string/1.4.4/lz-string.min.js
// @require      https://cdn.jsdelivr.net/npm/pako@1.0.11/dist/pako.min.js
// @grant        none
// ==/UserScript==

function mapModsByIndex(sequence, mods) {
    return mods.flatMap(mod => {
        // Find all indices of the amino acid in the sequence
        const indices = [];
        let idx = sequence.indexOf(mod.aminoAcid);
        while (idx !== -1) {
            indices.push(idx + 1); // Adding 1 to convert to 1-based index
            idx = sequence.indexOf(mod.aminoAcid, idx + 1);
        }

        // Return an array of objects for each index
        return indices.map(index => ({
            index: index,
            modMass: mod.modMass,
            aminoAcid: mod.aminoAcid
        }));
    }).filter(mod => mod !== null);
}



function mergeArrays(arr1, arr2) {
    return [...arr1, ...arr2];
}

function addModifications(sequence, varMods, ntermMod, ctermMod) {
    varMods.sort((a, b) => a.index - b.index);

    let addedLength = 0;
    varMods.forEach(mod => {
        let pos = mod.index + addedLength; // Adjust index for previously added characters
        if (pos <= sequence.length) {
            sequence = sequence.substring(0, pos) + "(" + mod.modMass + ")" + sequence.substring(pos);
            addedLength += mod.modMass.toString().length + 2; // Increment addedLength
        }
    });

    if (ntermMod !== 0.0) {
        sequence = "[" + ntermMod + "]" + sequence;
    }

    if (ctermMod !== 0.0) {
        sequence += "[" + ctermMod + "]";
    }

    return sequence;
}

function formatPeaks(peaksArray) {
    return peaksArray.map(function(peak) {
        const mzRounded = peak[0].toFixed(4);
        const intensityRounded = peak[1].toFixed(2);
        return mzRounded + ':' + intensityRounded;
    }).join(';');
}

function createStreamlitLink(baseURL, sequence, spectra) {
    let queryParams = new URLSearchParams();
    queryParams.set('sequence', sequence || 'DEFAULT_PEPTIDE');
    queryParams.set('spectra', spectra || 'DEFAULT_SPECTRA');
    return baseURL + '?' + queryParams.toString()
}

function addLinkToPage(url) {
    var linkDiv = document.createElement('div');
    linkDiv.style.cssText = `
            text-align: center;
            margin: 20px 0;
            padding: 10px 0;
            background-color: #f8f9fa;
            border: 1px solid #007bff;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0, 123, 255, 0.3);
        `;

    var link = document.createElement('a');
    link.href = url;
    link.textContent = 'Open Streamlit App';
    link.style.cssText = `
            color: #007bff;
            text-decoration: none;
            font-weight: bold;
            font-size: 16px;
            padding: 8px 15px;
        `;

    // Adding hover effect
    link.onmouseover = function() {
        this.style.backgroundColor = '#007bff';
        this.style.color = '#fff';
    };
    link.onmouseout = function() {
        this.style.backgroundColor = '';
        this.style.color = '#007bff';
    };

    linkDiv.appendChild(link);
    document.body.appendChild(linkDiv);
}

function floatToHex(f) {
    const buffer = new ArrayBuffer(4);
    new Float32Array(buffer)[0] = f;
    return [...new Uint8Array(buffer)].reverse().map(b => b.toString(16).padStart(2, '0')).join('');
}

function hexToFloat(hex) {
    const uint = new Uint32Array(new Uint8Array(hex.match(/.{1,2}/g).map(byte => parseInt(byte, 16))).buffer)[0];
    return new Float32Array([uint])[0];
}

function encodeLeadingZero(lz) {
    if (lz >= 0 && lz < 16) {
        return lz.toString(16);
    } else {
        throw new Error("Value must be between 0 and 15 (inclusive)");
    }
}

function reverseHexString(hexString) {
    return hexString.match(/.{2}/g).reverse().join('');
}

function hexDelta(a, b) {
    const diff = parseInt(a, 16) - parseInt(b, 16);
    return (diff >>> 0).toString(16).padStart(8, '0');
}

function countLeadingZeros(str) {
    return str.match(/^0*/)[0].length;
}

function encodeMzs(mzs) {
    const mzsHex = mzs.map(floatToHex);
    const initialHexValue = mzsHex[0];
    const initialHexValueZeros = countLeadingZeros(initialHexValue);
    const mzsHexDeltas = mzsHex.slice(1).map((hex, i) => hexDelta(hex, mzsHex[i]));
    const leadingZeros = mzsHexDeltas.map(countLeadingZeros);
    const hexDeltaStr = initialHexValue.slice(initialHexValueZeros) +
                        mzsHexDeltas.map((hex, i) => hex.slice(leadingZeros[i])).join('');
    const leadingZeroStr = encodeLeadingZero(initialHexValueZeros) +
                           leadingZeros.map(encodeLeadingZero).join('');
    return hexDeltaStr + leadingZeroStr.split('').reverse().join('');
}


function encodeIntensities(intensities) {
    return intensities.map(floatToHex).join('');
}

function splitPeaks(peaksArray) {
    const mzArray = peaksArray.map(peak => parseFloat(peak[0]));
    const intensityArray = peaksArray.map(peak => parseFloat(peak[1]));
    return [mzArray, intensityArray];
}

function CompressUrlLzstring(mzs, intensities) {
    let mz_str = encodeMzs(mzs);
    let int_str = encodeIntensities(intensities);
    const combinedData = JSON.stringify([mz_str, int_str]);
    return LZString.compressToEncodedURIComponent(combinedData);
}


const [mzArray, intensityArray] = splitPeaks(peaks);
console.log(mzArray.length);
const compressedData = CompressUrlLzstring(mzArray, intensityArray);
const allMods = mergeArrays(mapModsByIndex(mysequence, staticMods), varMods)
var modifiedSequence = addModifications(mysequence, allMods, ntermMod, ctermMod);

var localURL = 'http://localhost:8501/'
var liveURL = 'https://spectrum-viewer.streamlit.app/'
var streamlitURL = createStreamlitLink(liveURL, modifiedSequence, compressedData);
console.log(streamlitURL);
addLinkToPage(streamlitURL);
