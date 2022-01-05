// Sources:
// https://stackoverflow.com/questions/39479090/read-n-lines-of-a-big-text-file
// https://flask.palletsprojects.com/en/1.1.x/patterns/jquery/
// https://stackoverflow.com/questions/6831918/node-js-read-a-text-file-into-an-array-each-line-an-item-in-the-array/12299566
// https://stackoverflow.com/questions/6831918/node-js-read-a-text-file-into-an-array-each-line-an-item-in-the-array/12299566
// https://www.freecodecamp.org/news/javascript-from-callbacks-to-async-await-1cc090ddad99/
// https://developer.mozilla.org/de/docs/Web/JavaScript/Reference/Statements/async_function
// https://simon-schraeder.de/posts/filereader-async/
// Code copied and modified from js.js

class TextReader {
    CHUNK_SIZE = 8192000; // I FOUND THIS TO BE BEST FOR MY NEEDS, CAN BE ADJUSTED
    position = 0;
    length = 0;

    byteBuffer = new Uint8Array(0);

    lines = [];
    lineCount = 0;
    lineIndexTracker = 0;

    fileReader = new FileReader();
    textDecoder = new TextDecoder(`utf-8`);

    get allCachedLinesAreDispatched() {
        return !(this.lineIndexTracker < this.lineCount);
    }

    get blobIsReadInFull() {
        return !(this.position < this.length);
    }

    get bufferIsEmpty() {
        return this.byteBuffer.length === 0;
    }

    get endOfStream() {
        return this.blobIsReadInFull && this.allCachedLinesAreDispatched && this.bufferIsEmpty;
    }

    constructor(blob) {
        this.blob = blob;
        this.length = blob.size;
    }

    blob2arrayBuffer(blob) {
        return new Promise((resolve, reject) => {
            this.fileReader.onerror = reject;
            this.fileReader.onload = () => {
                resolve(this.fileReader.result);
            };

            this.fileReader.readAsArrayBuffer(blob);
        });
    }

    read(offset, count) {
        return new Promise(async (resolve, reject) => {
            if (!Number.isInteger(offset) || !Number.isInteger(count) || count < 1 || offset < 0 || offset > this.length - 1) {
                resolve(new ArrayBuffer(0));
                return
            }

            let endIndex = offset + count;

            if (endIndex > this.length) endIndex = this.length;

            let blobSlice = this.blob.slice(offset, endIndex);

            resolve(await this.blob2arrayBuffer(blobSlice));
        });
    }

    readLine() {
        return new Promise(async (resolve, reject) => {

            if (!this.allCachedLinesAreDispatched) {
                resolve(this.lines[this.lineIndexTracker++] + `\n`);
                return;
            }

            while (!this.blobIsReadInFull) {
                let arrayBuffer = await this.read(this.position, this.CHUNK_SIZE);
                this.position += arrayBuffer.byteLength;

                let tempByteBuffer = new Uint8Array(this.byteBuffer.length + arrayBuffer.byteLength);
                tempByteBuffer.set(this.byteBuffer);
                tempByteBuffer.set(new Uint8Array(arrayBuffer), this.byteBuffer.length);

                this.byteBuffer = tempByteBuffer;

                let lastIndexOfLineFeedCharacter = this.byteBuffer.lastIndexOf(10); // LINE FEED CHARACTER (\n) IS ONE BYTE LONG IN UTF-8 AND IS 10 IN ITS DECIMAL FORM

                if (lastIndexOfLineFeedCharacter > -1) {
                    let lines = this.textDecoder.decode(this.byteBuffer).split(`\n`);
                    this.byteBuffer = this.byteBuffer.slice(lastIndexOfLineFeedCharacter + 1);

                    let firstLine = lines[0];

                    this.lines = lines.slice(1, lines.length - 1);
                    this.lineCount = this.lines.length;
                    this.lineIndexTracker = 0;

                    resolve(firstLine + `\n`);
                    return;
                }
            }

            if (!this.bufferIsEmpty) {
                let line = this.textDecoder.decode(this.byteBuffer);
                this.byteBuffer = new Uint8Array(0);
                resolve(line);
                return;
            }

            resolve(null);
        });
    }
}

async function extract(){
    let file = document.getElementById("infile").files[0];
    let textReader = new TextReader(file);
    var reads = [];
    var lineno = 1;
    while (!textReader.endOfStream) {
        let line = await textReader.readLine();
        line = line.replace(/(\r\n|\n|\r)/gm, "");
        // only using the sequence reads

        //fasta file: taking all lines execpt the ones with '>'
        if ((line.charAt(0) !== '>') && (line != "")){
             reads.push(line);
        }else{
                reads.push('>');
            }

        // stopping after MAX(1Mio) lines
        if (reads.length == 1000000){
            break;
        }
        lineno++;
    }

    return reads
}

async function asyncCall(ext) {
  // Reading File
  const result = await extract();
  return result;
}

// Source:
// https://stackoverflow.com/questions/53694709/passing-javascript-array-in-python-flask

$(document).ready(function () {
    $("#remover").on("click", async function() {
        // prevent default send
        event.preventDefault();
        let name = document.getElementById("fname_rem").value;
        let valid = added.includes(name);

        if (!valid){
            alert('Error - Removing Filter: Please enter one of the deletable Filters!');
            return;
        }
        document.getElementById("content").style.display = "none";
        document.getElementById("loading-display").style.display = "block";
        var js_data = JSON.stringify(['REMOVE', name]);

        $.ajax({
        url: '/add_and_remove',
        type : 'post',
        contentType: 'application/json',
        dataType : 'json',
        data : js_data,
        success: function(){
        window.location.href = '/'
        }
     });
    });
});

$(document).ready(function () {
    $("#adder").on("click", async function() {
        // prevent default send
        event.preventDefault();

        // Testing if name is valid
        let name = document.getElementById("fname_add").value;
        var not_allowed = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5','IC6','IC7','IC8', 'None'].concat(added)
        let valid = not_allowed.includes(name);
        if ((valid) || (name.length < 3)){
            alert('Error - Adding Filter: Name already used or not at least 3 letters long!');
            return;
        }

        // Getting lines from File, max is 250.000
        let file = document.getElementById("infile").files[0];
        if (!file) {
            alert('No file selected, please select a fasta or fna file');
            return;
        }
        fname = document.getElementById("infile").files[0].name;
        ext = fname.split('.').pop();
        if ((ext !== 'fasta') && (ext !== 'fna')){
            alert('Error - Adding Filter: Wrong file-type, please select a fasta/fna file');
            return;
        }
        var reads = await asyncCall();

        // Getting new SVM DATA
        var lines = document.getElementById("SVM_add").value.split("\n");
        // Splitting Up Values
        for (var i = 0; i < (lines.length); i++) {
            lines[i] = lines[i].split(",");
        }

        // Each Filter needs at least 1 Vector
        if (lines.length < row + 1){
            alert('Error - Adding SVM: Number of Rows must be at least' + (row + 1) + ' !');
            return;
        }
        // Error Catching or Invalid Format Testing
        // First for Loop checks number of elements
        var allowed = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5','IC6','IC7','IC8', 'None'].concat(added)
        allowed.push(name);

        var labels = []
        for (var i = 0; i < (lines.length); i++) {
            if (lines[i].length != col + 1){
                alert('Error - Adding SVM: Row ' + (i+1) + ' has an invalid number of elements!');
                return;
            }
            //Last element must be label of existing filter
            var label = lines[i][lines[i].length - 1];
            var valid2 = allowed.includes(label);
            if (!valid2){
                alert('Error - Adding SVM: Row ' + (i+1) + ' has a Label of a not existing filter: ' + label);
                return;
            }
            labels.push(label);

            // Second for loop tests Values
             for (var j = 1; j < (lines[i].length - 1); j++){
                //Checking if [1:-1] are float
                val = parseFloat(lines[i][j]);
                if(isNaN(val)){
                    alert('Error - Adding SVM: Row ' + (i+1) + ', Col ' + j + ' must be a float!');
                    return;
                }
                // Checking Scores: Must be between 0 and 1
                if((val < 0)||(val>1)){
                    alert('Error - Adding SVM: Row ' + (i+1) + ', Col ' + j + ' must be between 0 and 1!');
                    return;
                }else{
                    lines[i][j] = val;
                }
             }
        }
        // Checking if each filter has at least one Vector
        var checker = allowed.every(elem => labels.includes(elem));
        if (!checker){
            alert('Error - Editing SVM: Each Filter must have at least 1 vector!');
            return;
        }


        document.getElementById("content").style.display = "none";
        document.getElementById("loading-display").style.display = "block";
        var js_data = JSON.stringify(['ADD', name, lines, reads]);


        $.ajax({
        url: '/add_and_remove',
        type : 'post',
        contentType: 'application/json',
        dataType : 'json',
        data : js_data,
        success: function(){
        window.location.href = '/'
        }
     });
    });
});


$(document).ready(function () {
    $("#SVM").on("click", async function() {
        event.preventDefault();
        var lines = document.getElementById("SVM_TA").value.split("\n");
        // Splitting Up Values
        for (var i = 0; i < (lines.length); i++) {
            lines[i] = lines[i].split(",");
        }

        // Every Filter needs at least one Vector
        if (lines.length < row_min){
            alert('Error - Editing SVM: Number of Rows must be at least ' + row_min + ' !');
            return;
        }
        // Error Catching or Invalid Format Testing
        // First for Loop checks number of elements
        var allowed = ['IC1', 'IC2', 'IC3', 'IC4', 'IC5','IC6','IC7','IC8', 'None'].concat(added)
        var labels = []
        for (var i = 0; i < (lines.length); i++) {
            if (lines[i].length != col){
                alert('Error - Editing SVM: Row ' + (i+1) + ' has an invalid number of elements!');
                return;
            }
            //Last element must be label of existing filter
            var label = lines[i][lines[i].length - 1];
            var valid = allowed.includes(label);
            if (!valid){
                alert('Error - Editing SVM: Row ' + (i+1) + ' has a Label of a not existing filter: ' + label);
                return;
            }
            labels.push(label);

            // Second for loop tests Values
             for (var j = 1; j < (lines[i].length - 1); j++){
                //Checking if [1:-1] are float
                val = parseFloat(lines[i][j]);
                if(isNaN(val)){
                    alert('Error - Editing SVM: Row ' + (i+1) + ', Col ' + j + ' must be a float!');
                    return;
                }
                // Checking Scores: Must be between 0 and 1
                if((val < 0)||(val>1)){
                    alert('Error - Editing SVM: Row ' + (i+1) + ', Col ' + j + ' must be between 0 and 1!');
                    return;
                }else{
                    lines[i][j] = val;
                }
             }
        }
        // Checking if at least 1 vector per filter
        // Source:
        // https://stackoverflow.com/questions/53606337/check-if-array-contains-all-elements-of-another-array
        var checker = allowed.every(elem => labels.includes(elem));
        if (!checker){
            alert('Error - Editing SVM: Each Filter must have at least 1 vector!');
            return;
        }

        document.getElementById("content").style.display = "none";
        document.getElementById("loading-display").style.display = "block";
        var js_data = JSON.stringify(['SVM', lines]);

        $.ajax({
        url: '/add_and_remove',
        type : 'post',
        contentType: 'application/json',
        dataType : 'json',
        data : js_data,
        success: function(){
        window.location.href = '/'
        }
     });
    });
});

$(document).ready(function () {
    $("#rem_oxa").on("click", async function() {
        // prevent default send
        event.preventDefault();
        let name = document.getElementById("oxa_id_rem").value;
        let valid = oxa_ids.includes(name)
        if (!valid){
            alert('Error - Removing OXA: Please enter one of the deletable OXAs!');
            return;
        }
        document.getElementById("content").style.display = "none";
        document.getElementById("loading-display").style.display = "block";
        var js_data = JSON.stringify(['REMOVE_OXA', name]);

        $.ajax({
        url: '/add_and_remove',
        type : 'post',
        contentType: 'application/json',
        dataType : 'json',
        data : js_data,
        success: function(){
        window.location.href = '/'
        }
     });
    });
});

$(document).ready(function () {
    $("#add_oxa").on("click", async function() {
        // prevent default send
        event.preventDefault();
        let name = document.getElementById("oxa_id").value;
        let valid = oxa_ids.includes(name)
        if ((valid) || (name.length < 3)){
            alert('Error - ADD OXA: Name already in use or contains less than 3 Chars!');
            return;
        }

        let file = document.getElementById("oxa_file").files[0];
        if (!file) {
            alert('OXA-Adder: No file selected, please select a fasta file');
            return;
        }
        filen = document.getElementById("oxa_file").files[0].name;
        ext = filen.split('.').pop();

        if (ext !== 'fasta'){
            alert('OXA-Adder: Wrong file-type, please select a fasta-file!');
            return;
        }

        document.getElementById("content").style.display = "none";
        document.getElementById("loading-display").style.display = "block";
        var reads = await asyncCall_oxa();
        var js_data = JSON.stringify(['ADD_OXA', name, reads]);

        $.ajax({
        url: '/add_and_remove',
        type : 'post',
        contentType: 'application/json',
        dataType : 'json',
        data : js_data,
        success: function(){
        window.location.href = '/'
        }
     });
    });
});

async function get_oxa(){
    let file = document.getElementById("oxa_file").files[0];
    let textReader = new TextReader(file);
    var reads = [];
    var lineno = 1;
    while (!textReader.endOfStream) {
        let line = await textReader.readLine();
        line = line.replace(/(\r\n|\n|\r)/gm, "");
        // only using the sequence reads

        //fasta file: taking all lines execpt the ones with '>'
        if ((line.charAt(0) !== '>') && (line != "")){
             reads.push(line);
        }else{
                reads.push('>');
            }

        // stopping after MAX(50000) lines
        if (reads.length == 50000){
            break;
        }
        lineno++;
    }

    return reads
}

async function asyncCall_oxa(ext) {
  // Reading File
  const result = await get_oxa();
  return result;
}
