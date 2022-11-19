// Sources:
// https://stackoverflow.com/questions/39479090/read-n-lines-of-a-big-text-file
// https://flask.palletsprojects.com/en/1.1.x/patterns/jquery/
// https://stackoverflow.com/questions/6831918/node-js-read-a-text-file-into-an-array-each-line-an-item-in-the-array/12299566
// https://stackoverflow.com/questions/6831918/node-js-read-a-text-file-into-an-array-each-line-an-item-in-the-array/12299566
// https://www.freecodecamp.org/news/javascript-from-callbacks-to-async-await-1cc090ddad99/
// https://developer.mozilla.org/de/docs/Web/JavaScript/Reference/Statements/async_function
// https://simon-schraeder.de/posts/filereader-async/
class TextReader {
    // https://stackoverflow.com/a/55377748/9100798
    CHUNK_SIZE = 8192000;
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

async function extract(Max_reads, checks, ext){
    let file = document.getElementById("infile").files[0];
    let textReader = new TextReader(file);
    var name = document.getElementById('infile').files[0].name;
    var reads = [];
    var lineno = 1;
    var max = Max_reads;
    while (!textReader.endOfStream) {
        let line = await textReader.readLine();
        line = line.replace(/(\r\n|\n|\r)/gm, "");
        // only using the sequence reads
        if ((ext === 'fq')||(ext === 'fastq')){
            if (((lineno)%2==0)&&((lineno)%4!=0)){
                reads.push(line);
            }
        }else{
            //fasta file: taking all lines
            if (line.charAt(0) !== '>'){
                reads.push(line);
            }
            else{
                reads.push('>');
            }
        }

        // stopping after n lines
        if (reads.length == max){
            break;
        }

        lineno++;
    }
    reads.push(name);
    result = reads.concat(checks)
    return result
}

async function asyncCall(ext) {
  document.getElementById("extracter").style.display = "block";
  var max_reads = document.getElementById("reads_max").value;
  var checks = [];

  // saving checkbox info so the Form can be hidden
  checks.push(document.getElementById("quick").checked);
  checks.push(document.getElementById("IC1").checked);
  checks.push(document.getElementById("IC2").checked);
  checks.push(document.getElementById("IC3").checked);
  checks.push(document.getElementById("IC4").checked);
  checks.push(document.getElementById("IC5").checked);
  checks.push(document.getElementById("IC6").checked);
  checks.push(document.getElementById("IC7").checked);
  checks.push(document.getElementById("IC8").checked);
  checks.push(document.getElementById("added").checked);
  checks.push(document.getElementById("OXA").checked);
  // Deactivating Checkboxes etc while extracting reads
  document.getElementById("opt").style.display = "none";

  // Complete fileupload (Max 100.000 lines) if fasta file
  if ((ext == 'fasta') || (ext == 'fna')){
    max_reads = 100000;
  }
  if (((ext == 'fq') || (ext == 'fastq')) && (document.getElementById("OXA").checked)){
    max_reads = 250000;
  }
  // Assigning
  const result = await extract(max_reads, checks, ext);

  return result;
}

async function asyncCallspec(ext) {
  document.getElementById("extracter").style.display = "block";
  var max_reads = document.getElementById("reads_max").value;
  var checks = [];

  // saving checkbox info so the Form can be hidden
  checks.push(document.getElementById("quick").checked);
  checks.push(document.getElementById("OXA").checked);
  checks.push(document.getElementById("metagenome").checked);
  //checks.push(document.getElementById("added").checked);
  // Deactivating Checkboxes etc while extracting reads
  document.getElementById("opt").style.display = "none";

  // Complete fileupload (Max 100.000 lines) if fasta file
  if ((ext == 'fasta') || (ext == 'fna')){
    max_reads = 5000000;
  }
  //if (((ext == 'fq') || (ext == 'fastq')) && (document.getElementById("OXA").checked)){
  //  max_reads = 250000;
  //}

  // Assigning
  const result = await extract(max_reads, checks, ext);

  return result;
}

// Source:
// https://stackoverflow.com/questions/53694709/passing-javascript-array-in-python-flask
$(document).ready(function () {
    $("#submit").on("click", async function() {
        // prevent default send
        event.preventDefault();

        let file = document.getElementById("infile").files[0];
        if (!file) {
            alert('No file selected, please select a .fq file that contains sequence reads');
            return;
        }
        name = document.getElementById('infile').files[0].name;
        ext = name.split('.').pop();

        if ((ext !== 'fq') && (ext !== 'fasta') && (ext !== 'fna') && (ext !== 'fastq')){
            alert('Wrong file-type, please select a FASTQ or FASTA/FNA file');
            return;
        }

        var number = document.getElementById("reads_max").value;

        // Converting Number field to String, then checking if
        // only numbers are in string (to prevent entering '.' or '+'
        var not_int = !(/^\d+$/.test(number.toString()));

        if ((not_int) || (number < 500) || (number > 1000000)){
            alert('Error: Number of reads must be between 500 and 100.000 and also be a Integer!');
            return;
        }

        // Getting Reads
        var js_data = JSON.stringify(await asyncCall(ext));

        if (js_data == null){
            alert('Error: This Tool does not support your Browser, please use a modern Browser.');
            return;
        }

        $.ajax({
            url: '/ic',
            type : 'post',
            contentType: 'application/json',
            dataType : 'json',
            data : js_data,
            success: function(){
          //  document.getElementById("content").style.display = "none";
        //    document.getElementById("loading-display").style.display = "block";
            window.location.href = '/assign'

         },
            error: function() {
              //  document.getElementById("opt").style.display = "block";
              //  document.getElementById("extracter").style.display = "none";
              //  document.getElementById("loading-display").style.display = "none";
                alert("Your Browser does not support this Tool. Please use a valid Browser");
            }
        });
    });
});



function myFunc(literature){
    for (let i = 0; i < literature[0].length; i++) {
        const test1 = document.createElement("li");
        const test2 = document.createElement("p");
        test2.style.lineHeight = "75%";
        const test3 = document.createElement("a");
        test3.setAttribute('href',literature[0][i]);
        test3.innerText = literature[1][i];
        const test4 = document.createElement("a");
        test4.className = "btn";
        test4.setAttribute("data-bs-toggle", "collapse");
        test4.setAttribute("data-bs-target", "#" + literature[5][i]);
        const test5 = document.createElement("small");
        test5.innerText = "[View Details]";
        const test6_1 = document.createElement("br");
        const test6_2 = document.createElement("br");
        const test6_3 = document.createElement("br");
        //const test6_4 = document.createElement("br");
        //const test6_5 = document.createElement("br");
        const test7 = document.createElement("small");
        const test8 = document.createElement("small");
        const test9 = document.createElement("font");
        test9.setAttribute('color',"4d8055");
        const test6 = document.createElement("br");
        const test10 = document.createElement("p");
        test10.className = "collapse";
        test10.setAttribute('id', literature[5][i]);
        const test11 = document.createElement("b");
        test11.innerText = "Abstract:"

        test10.appendChild(test11);
        test10.appendChild(test6);
        var newContent = document.createTextNode(literature[2][i]);
        test10.appendChild(newContent);

        var newContent = document.createTextNode(literature[4][i]);
        test9.appendChild(newContent);
        test8.appendChild(test9);
        var newContent = document.createTextNode(literature[3][i]);
        test7.appendChild(newContent);
        test4.appendChild(test5);

        test2.appendChild(test3);
        test2.appendChild(test4);
        test2.append(test6_1);
        test2.appendChild(test7);
        test2.append(test6_2);
        test2.append(test6_3);
      //  test2.append(test6_4);
        //test2.append(test6_5);
        test2.appendChild(test8);

        test1.appendChild(test2);
        test1.appendChild(test10);

        document.getElementById("literature").appendChild(test1);
      }
      return literature
    }


$(document).ready(function(){
  $("#close_popup").on("click", async function() {
    event.preventDefault();
    document.getElementById("popup-1").classList.toggle("active");
  });
});


$(document).ready(function(){
  $("#tab-3").on("click", async function() {
    event.preventDefault();
    if (document.getElementById("popup-2").classList == "popup-2 active") {
    document.getElementById("popup-2").classList.remove("active");
}
  });
});


$(document).ready(function(){
  $("#Display_options").on("click", async function() {
    event.preventDefault();
    document.getElementById("popup-2").classList.toggle("active");
  });
});



$(document).ready(function(){
        $("#infile").change(function(){
          name = document.getElementById('infile').files[0].name;
          if (document.getElementById("popup-1").classList == "popup") {
          document.getElementById("popup-1").classList.add("active");
    }
          ext = name.split('.').pop();
          if ((ext == 'fq') || (ext == 'fastq')){
              y = document.getElementById("AspecTinput");
              y.style.display = "block";
              }
          else {
              y = document.getElementById("AspecTinput");
              y.style.display = "none";
          }
          if ((ext == 'fasta') || (ext == 'fna')){
              y = document.getElementById("AspecTinput-2");
              y.style.display = "block";
              }
          else {
              y = document.getElementById("AspecTinput-2");
              y.style.display = "none";
          }
        });
    });


$(document).ready(function(literature) {
    $("#apply").on("click", async function() {
      // prevent default send
      event.preventDefault();
      document.getElementById("popup-2").classList.toggle("active");
    //  const clear_literature = document.getElementById("literature");
    //  clear_literature.innerHTML = '';

      var js_data1 = document.getElementById("literature_max").value;
      var js_data2 = document.getElementById("id_sort").value;
      const js_data = JSON.stringify([js_data1, js_data2]);

      $.ajax({
          url: '/resultsspec',
          type : 'post',
          contentType: 'application/json; charset=utf-8',
          dataType : 'json',
          data : js_data,
          success: function(data){
            const clear_literature = document.getElementById("literature");
            clear_literature.innerHTML = '';
            myFunc(data);

         },
          error: function() {
            const clear_literature = document.getElementById("literature");
            clear_literature.innerHTML = '';
            myFunc(data);
            }
      });
      });
    });


$(document).ready(function () {
    $("#submitspec").on("click", async function() {
        // prevent default send
        event.preventDefault();

        let file = document.getElementById("infile").files[0];
        if (!file) {
            alert('No file selected, please select a .fq file that contains sequence reads');
            return;
        }
        name = document.getElementById('infile').files[0].name;
        ext = name.split('.').pop();

        if ((ext !== 'fq') && (ext !== 'fasta') && (ext !== 'fna') && (ext !== 'fastq')){
            alert('Wrong file-type, please select a FASTQ or FASTA/FNA file');
            return;
        }

        var number = document.getElementById("reads_max").value;

        // Converting Number field to String, then checking if
        // only numbers are in string (to prevent entering '.' or '+'
        var not_int = !(/^\d+$/.test(number.toString()));

        if ((not_int) || (number < 5000) || (number > 10000000)){
            alert('Error: Number of reads must be between 5000 and 10.000.000 and also be a Integer!');
            return;
        }

        // Getting Reads
      //  var data = await asyncCallspec(ext)
      //  var js_data="[";
      //    for(var indx=0;indx<data.length-1;indx++){
      //      js_data+=JSON.stringify(data[indx],null,4)+",";
      //    }
      //    js_data+=JSON.stringify(data[data.length-1],null,4)+"]";
        var js_data = JSON.stringify(await asyncCallspec(ext));

        if (js_data == null){
            alert('Error: This Tool does not support your Browser, please use a modern Browser.');
            return;
        }


        $.ajax({
            url: '/species',
            type : 'post',
            contentType: 'application/json',
            dataType : 'json',
            data : js_data,
            success: function(){
            //document.getElementById("content").style.display = "none";
            //document.getElementById("loading-display").style.display = "block";
            window.location.href = '/assignspec';

         },
            error: function() {
                document.getElementById("opt").style.display = "block";
                document.getElementById("extracter").style.display = "none";
                document.getElementById("loading-display").style.display = "none";
                alert("Your Browser does not support this Tool. Please use a valid Browser");
            }
        });

    });
});
