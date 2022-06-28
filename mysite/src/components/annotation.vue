<template>
  <div id="annotation" style="color: rgb(33, 56, 50) !important">
    <div style="display: flex; justify-content: center">
      <p><strong>Annotation</strong>&nbsp;by SELINA</p>
    </div>
    <el-row>
      <el-col :span="24" style="padding=10px">
        <p style="font-weight: bold">
          <strong>
            <a href="https://selina.readthedocs.io/en/latest/index.html"
              >SELINA</a
            ></strong
          >
          is a deep learning-based framework for single cell assignment with
          multiple references.  
              <el-tooltip class="item" effect="dark" content="Click for instructions" placement="top-start">

          <el-button
          id = 'help-button'
          @click="drawer = true"
          type="primary"
          style="margin-left: 16px"
          icon="el-icon-chat-line-round"
          class="custom-icon"
          circle
          ></el-button
        >
        </el-tooltip>
        </p>
      </el-col>
    </el-row>
    <el-card class="anno">
      <el-row>
        <el-divider content-position="center"
          ><span style="color: rgb(26, 96, 67); font-weight: bold">Input</span>
        </el-divider>
      </el-row>

      <el-row :gutter="15" type="flex" style="align-items: center">
        <el-col :span="4" style="color: rgb(33, 56, 50)">
          <p>Count matrix file:</p>
        </el-col>
        <el-col :span="8">
          <form class="file">
            <input
              type="file"
              class="file-btn"
              value=""
              title=" "
              id="file"
              accept=".gz"
              @change="getFile_main($event)"
            />
          </form>
        </el-col>

        <el-col
          v-if="mode_value == 'single'"
          :span="4"
          style="color: rgb(33, 56, 50)"
        >
          <p>Matrix format:</p>
        </el-col>
        <el-col v-if="mode_value == 'single'" :span="4">
          <el-select v-model="type_value" placeholder="h5/mtx/plain">
            <el-option
              v-for="item in file_type"
              :key="item.value"
              :label="item.label"
              :value="item.value"
            >
            </el-option>
          </el-select>
        </el-col>
      </el-row>
      <el-row :gutter="15" type="flex" style="align-items: center">
        <el-col :span="4" style="color: rgb(33, 56, 50)">
          <p>Expression level:</p>
        </el-col>
        <el-col :span="4">
          <el-select v-model="mode_value" placeholder="Single/Cluster">
            <el-option
              v-for="item in mode"
              :key="item.value"
              :label="item.label"
              :value="item.value"
            >
            </el-option>
          </el-select>
        </el-col>
      </el-row>

      <el-row v-if="type_value == 'mtx'" :gutter="15">
        <el-col
          :span="4"
          style="color: rgb(33, 56, 50)"
          v-if="type_value == 'mtx'"
        >
          <p>Barcode:</p>
        </el-col>
        <el-col :span="8" v-if="type_value == 'mtx'">
          <form>
            <input
              type="file"
              value=""
              id="file"
              accept=".gz"
              @change="getFile_bar($event)"
              v-if="type_value == 'mtx'"
            />
          </form>
        </el-col>
        <el-col
          :span="4"
          style="color: rgb(33, 56, 50)"
          v-if="type_value == 'mtx'"
        >
          <p>Feature:</p>
        </el-col>
        <el-col :span="4" v-if="type_value == 'mtx'">
          <form>
            <input
              type="file"
              value=""
              id="file"
              accept=".gz"
              @change="getFile_gen($event)"
              v-if="type_value == 'mtx'"
            />
          </form>
        </el-col>
      </el-row>

      <el-row style="align-items: center">
        <el-divider
          ><span style="color: rgb(26, 96, 67); font-weight: bold"
            >Parameters</span
          >
        </el-divider>
        <el-row type="flex" style="align-items: center">
          <el-col :span="4" style="color: rgb(33, 56, 50)">
            <p>Model selection:</p>
          </el-col>
          <el-row style="align-items: center">
            <el-col :span="12">
              <el-select
                v-model="reference_value"
                placeholder="Normal/Disease"
                @change="getref"
              >
                <el-option
                  v-for="item in reference"
                  :key="item.value"
                  :label="item.label"
                  :value="item.value"
                >
                </el-option>
              </el-select>
            </el-col>
            <el-col :span="12">
              <el-select v-model="model_value" placeholder="Tissue">
                <el-option
                  v-for="item in model"
                  :key="item.value"
                  :label="item.label"
                  :value="item.value"
                  :disabled="isDisabled"
                >
                </el-option>
              </el-select>
            </el-col>
          </el-row>
          <el-col
            :span="4"
            style="color: rgb(33, 56, 50)"
            v-if="type_value == 'plain'"
          >
            <p>Separator:</p>
          </el-col>
          <el-col :span="4" v-if="type_value == 'plain'">
            <el-select v-model="separator_value" placeholder="Tab/Space/Comma">
              <el-option
                v-for="item in separator"
                :key="item.value"
                :label="item.label"
                :value="item.value"
              >
              </el-option>
            </el-select>
          </el-col>
        </el-row>
        <el-row type="flex" style="align-items: center">
          <el-col :span="4" style="color: rgb(33, 56, 50)">
            <p>Gene:</p>
          </el-col>
          <el-col :span="4">
            <el-select v-model="idtype_value" placeholder="Symbol/Ensembl">
              <el-option
                v-for="item in idtype"
                :key="item.value"
                :label="item.label"
                :value="item.value"
              >
              </el-option>
            </el-select>
          </el-col>
          <el-col :span="4" style="color: rgb(33, 56, 50)">
            <p>Assembly:</p>
          </el-col>
          <el-col :span="4">
            <el-select v-model="assembly_value" placeholder="GRCh38/37">
              <el-option
                v-for="item in assembly"
                :key="item.value"
                :label="item.label"
                :value="item.value"
              >
              </el-option>
            </el-select>
          </el-col>

          <el-col
            :span="4"
            style="color: rgb(33, 56, 50)"
            v-if="mode_value == 'single'"
          >
            <p>Clustering resolution:</p>
          </el-col>
          <el-col :span="4" v-if="mode_value == 'single'">
            <el-input v-model="res" placeholder="0.0~1.0"></el-input>
          </el-col>
        </el-row>
      </el-row>
      <el-row type="flex" style="align-items: center">
        <el-col
          :span="4"
          style="color: rgb(33, 56, 50)"
          v-if="mode_value == 'single'"
        >
          <p>PC number:</p>
        </el-col>
        <el-col :span="4" v-if="mode_value == 'single'">
          <el-input v-model="npc" placeholder="Dimensions after PCA"></el-input>
        </el-col>
        <el-col :span="4" style="color: rgb(33, 56, 50)">
          <p>Email:</p>
        </el-col>
        <el-col :span="4">
          <el-input
            v-model="email"
            placeholder="Email Address"
            prefix-icon="el-icon-message"
          >
          </el-input>
        </el-col>
      </el-row>

      <div style="display: flex; justify-content: center">
        <el-button type="primary" round @click="submitForm($event)"
          >Submmit</el-button
        >
      </div>
    </el-card>
    <!-- <el-row>
      <el-col :span="4" :offset="20">
        <el-button
          @click="drawer = true"
          type="primary"
          style="margin-left: 16px"
          circle
          >Help</el-button
        >
      </el-col>
    </el-row> -->

    <el-drawer :visible.sync="drawer">
      <h1 style="text-align: center">Instruction</h1>
      <div id="Format">
        <h3>Input</h3>
        <br />
        <p>
          <b>Count matrix file:</b> SELINA supports 3 input formats:
          <code>plain</code>, <code>h5</code> and <code>mtx</code>. You can find
          detail information in
          <a
            class="doclink"
            href="https://selina.readthedocs.io/en/latest/prepare.html#preprocess-of-query-data"
            target="_blank"
            >Preprocess-of-query-data</a
          >
          part.
        </p>
        <p>
          To guarantee your experience, we only allow for
          <code>.tar.gz</code> files for time saving. You need to compress the
          original file. Noted that the compressed file should have same prefix
          with the original file name (eg :
          <code>tar -zcvf Query.txt.tar.gz Query.txt</code>).
        </p>
        <br />
        <p>
          <b>Expression level:</b> It depends on the expression level of your
          matrix uploaded, choose from single-cell level or cluster level.
        </p>
        <br />
        <p>
          <b>Matrix format:</b> It depends on which format your input file is
          (only need when expression level is 'single').
        </p>
        <br />
        <h3>Parameters</h3>
        <br />
        <p>
          <b>Model Selection: </b> You need to choose data type from
          <code>normal tissue</code> or <code>disease tissue</code>, which
          depends on the condition your input is. For tissue input, we have
          trained MADA models for 35 kinds of tissues, five of which also have
          correspoding disease model(Intestine, Pancreas, Skin, Lung, Liver).
          You can choose tissue model based on your data source.
        </p>

        <br />
        <p>
          <b>Gene:</b> Type of gene name, <code>symbol</code> for gene symbol
          and <code> ensembl</code> for ensembl id.
        </p>
        <br />
        <p>
          <b>Assembly:</b> Genome version of the genes (GRCh38/hg38 or
          GRCh37/hg19).
        </p>
        <br />
        <p>
          <b>Separator:</b> It depends on which delimiters to separate text
          strings in your input (only need when input format is 'plain').
        </p>
        <br />
        <p>
          <b>Barcode:</b> Upload the barcode file (required for the input format
          of 'mtx').
        </p>
        <br />
        <p>
          <b>Feature:</b> Upload the feature file (required for the input format
          of 'mtx').
        </p>
        <br />

        <p><b>PC number:</b> Number of dimensions after PCA.</p>
        <br />
        <p>
          <b>Clustering resolution:</b> Resolution for the clustering process.
        </p>
        <br />

        <p>
          <b>Email Address:</b> We will send you annotation results after things
          done.
        </p>
        <br />

        <h3>Give it a try!</h3>
        <br />
        <p>
          You can download example file for testing, please select
          <code>plain</code> for format and <code>Pancreas</code> for tissue. We
          only kept only has 400 cells in example file to avoid network
          congestion. To ensure preprocess module works, we recommend you choose
          PC number smaller than <code>5</code>.

          <el-button id="dw" type="primary" round @click="download()"
            >Example File</el-button
          >
        </p>
        <br />
        <h3>More</h3>
        <br />

        <li>
          The data will be deleted after processing, please don't worry about
          your data security.
        </li>
        <li>
          If you want to annotate single-cell data from tissues we haven't
          collect, you can train your own model and play with detailed
          parameters settings using the local version of
          <a
            class="doclink"
            href="https://selina.readthedocs.io/en/latest/run.html"
            target="_blank"
            >SELINA</a
          >.
        </li>
        <li>
          Due to the limit of computational resource, it will take a while to
          get results for large dataset. If you can not find result email in
          your normal box, it is worth checking in your spam or junk mail
          section.
        </li>
      </div>
    </el-drawer>
  </div>
</template>

<script>
export default {
  name: "upload",
  data() {
    return {
      name: "UploadFile",
      checkStatus: "",
      mode_value: "",
      model_value: "",
      type_value: "",
      email: "",
      idtype_value: "",
      assembly_value: "",
      reference_value: "",
      separator_value: "",
      npc: "",
      res: "",
      drawer: false,
      labelPosition: "left",
      formLabelAlign: {
        name: "",
        region: "",
        type: "",
      },

      file_type: [
        { value: "h5", label: "h5" },
        { value: "mtx", label: "mtx" },
        { value: "plain", label: "plain" },
      ],
      mode: [
        { value: "single", label: "Single" },
        { value: "cluster", label: "Cluster" },
      ],
      idtype: [
        { value: "symbol", label: "Symbol" },
        { value: "ensembl", label: "Ensembl" },
      ],
      assembly: [
        { value: "GRCh38", label: "GRCh38" },
        { value: "GRCh37", label: "GRCh37" },
      ],
      reference: [
        { value: "normal", label: "Normal tissue" },
        { value: "disease", label: "Disease tissue" },
      ],
      separator: [
        { value: "tab", label: "Tab" },
        { value: "space", label: "Space" },
        { value: "comma", label: "Comma" },
      ],
      isDisabled: true,
      model: [
        { value: "Adrenal-Gland", label: "Adrenal-Gland" },
        { value: "Airway-Epithelium", label: "Airway-Epithelium" },
        { value: "Artery", label: "Artery" },
        { value: "Bladder", label: "Bladder" },
        { value: "Blood", label: "Blood" },
        { value: "Bone-Marrow", label: "Bone-Marrow" },
        { value: "Brain", label: "Brain" },
        { value: "Breast", label: "Breast" },
        { value: "Choroid", label: "Choroid" },
        { value: "Decidua", label: "Decidua" },
        { value: "Esophagus", label: "Esophagus" },
        { value: "Eye", label: "Eye" },
        { value: "Fallopian", label: "Fallopian" },
        { value: "Gall-Bladder", label: "Gall-Bladder" },
        { value: "Heart", label: "Heart" },
        { value: "Intestine", label: "Intestine" },
        { value: "Kidney", label: "Kidney" },
        { value: "Liver", label: "Liver" },
        { value: "Lung", label: "Lung" },
        { value: "Muscle", label: "Muscle" },
        { value: "Nose", label: "Nose" },
        { value: "Ovary", label: "Ovary" },
        { value: "Pancreas", label: "Pancreas" },
        { value: "Peritoneum", label: "Peritoneum" },
        { value: "Placenta", label: "Placenta" },
        { value: "Pleura", label: "Pleura" },
        { value: "Prostate", label: "Prostate" },
        { value: "Skin", label: "Skin" },
        { value: "Spleen", label: "Spleen" },
        { value: "Stomach", label: "Stomach" },
        { value: "Testis", label: "Testis" },
        { value: "Thyroid", label: "Thyroid" },
        { value: "Ureter", label: "Ureter" },
        { value: "Uterus", label: "Uterus" },
        { value: "Visceral-Adipose", label: "Visceral-Adipose" },
      ],
    };
  },

  methods: {
    getref(code) {
      if (code === "disease") {
        this.isDisabled = false;
        this.model = [
          { value: "Intestine", label: "Intestine" },
          { value: "Pancreas", label: "Pancreas" },
          { value: "Skin", label: "Skin" },
          { value: "Lung", label: "Lung" },
          { value: "Liver", label: "Liver" },
        ];
      } else if (code === "normal") {
        this.isDisabled = false;
        this.model = [
          { value: "Adrenal-Gland", label: "Adrenal-Gland" },
          { value: "Airway-Epithelium", label: "Airway-Epithelium" },
          { value: "Artery", label: "Artery" },
          { value: "Bladder", label: "Bladder" },
          { value: "Blood", label: "Blood" },
          { value: "Bone-Marrow", label: "Bone-Marrow" },
          { value: "Brain", label: "Brain" },
          { value: "Breast", label: "Breast" },
          { value: "Choroid", label: "Choroid" },
          { value: "Decidua", label: "Decidua" },
          { value: "Esophagus", label: "Esophagus" },
          { value: "Eye", label: "Eye" },
          { value: "Fallopian", label: "Fallopian" },
          { value: "Gall-Bladder", label: "Gall-Bladder" },
          { value: "Heart", label: "Heart" },
          { value: "Intestine", label: "Intestine" },
          { value: "Kidney", label: "Kidney" },
          { value: "Liver", label: "Liver" },
          { value: "Lung", label: "Lung" },
          { value: "Muscle", label: "Muscle" },
          { value: "Nose", label: "Nose" },
          { value: "Ovary", label: "Ovary" },
          { value: "Pancreas", label: "Pancreas" },
          { value: "Peritoneum", label: "Peritoneum" },
          { value: "Placenta", label: "Placenta" },
          { value: "Pleura", label: "Pleura" },
          { value: "Prostate", label: "Prostate" },
          { value: "Skin", label: "Skin" },
          { value: "Spleen", label: "Spleen" },
          { value: "Stomach", label: "Stomach" },
          { value: "Testis", label: "Testis" },
          { value: "Thyroid", label: "Thyroid" },
          { value: "Ureter", label: "Ureter" },
          { value: "Uterus", label: "Uterus" },
          { value: "Visceral-Adipose", label: "Visceral-Adipose" },
        ];
      } else {
        this.model = [
          { value: "Adrenal-Gland", label: "Adrenal-Gland" },
          { value: "Airway-Epithelium", label: "Airway-Epithelium" },
          { value: "Artery", label: "Artery" },
          { value: "Bladder", label: "Bladder" },
          { value: "Blood", label: "Blood" },
          { value: "Bone-Marrow", label: "Bone-Marrow" },
          { value: "Brain", label: "Brain" },
          { value: "Breast", label: "Breast" },
          { value: "Choroid", label: "Choroid" },
          { value: "Decidua", label: "Decidua" },
          { value: "Esophagus", label: "Esophagus" },
          { value: "Eye", label: "Eye" },
          { value: "Fallopian", label: "Fallopian" },
          { value: "Gall-Bladder", label: "Gall-Bladder" },
          { value: "Heart", label: "Heart" },
          { value: "Intestine", label: "Intestine" },
          { value: "Kidney", label: "Kidney" },
          { value: "Liver", label: "Liver" },
          { value: "Lung", label: "Lung" },
          { value: "Muscle", label: "Muscle" },
          { value: "Nose", label: "Nose" },
          { value: "Ovary", label: "Ovary" },
          { value: "Pancreas", label: "Pancreas" },
          { value: "Peritoneum", label: "Peritoneum" },
          { value: "Placenta", label: "Placenta" },
          { value: "Pleura", label: "Pleura" },
          { value: "Prostate", label: "Prostate" },
          { value: "Skin", label: "Skin" },
          { value: "Spleen", label: "Spleen" },
          { value: "Stomach", label: "Stomach" },
          { value: "Testis", label: "Testis" },
          { value: "Thyroid", label: "Thyroid" },
          { value: "Ureter", label: "Ureter" },
          { value: "Uterus", label: "Uterus" },
          { value: "Visceral-Adipose", label: "Visceral-Adipose" },
        ];
      }
    },

    download() {
      let link = document.createElement("a");
      link.setAttribute("download", "");
      link.href = "./backend/download/example.txt.tar.gz"; // 你本地资源文件的存放地址
      console.log("href:", link.href);
      link.click();
    },

    getFile_main(event) {
      (this.file_main = event.target.files[0]), console.log(this.file_main);
    },
    getFile_bar(event) {
      (this.file_bar = event.target.files[0]), console.log(this.file_bar);
    },
    getFile_gen(event) {
      (this.file_gen = event.target.files[0]), console.log(this.file_gen);
    },
    submitForm(event) {
      // alert('hi')
      // event.preventDefault();
      let formData = new FormData();
      formData.append("file_main", this.file_main);
      formData.append("file_bar", this.file_bar);
      formData.append("file_gen", this.file_gen);
      formData.append("mode", this.mode_value);
      formData.append("model", this.model_value);
      formData.append("email", this.email);
      formData.append("format", this.type_value);
      formData.append("idtype", this.idtype_value);
      formData.append("assembly", this.assembly_value);
      formData.append("reference", this.reference_value);
      formData.append("separator", this.separator_value);
      formData.append("npc", this.npc);
      formData.append("res", this.res);
      console.log(formData.get("reference"));
      console.log(formData.get("separator"));
      let config = {
        headers: {
          "Content-Type": "multipart/form-data",
        },
      };

      // 创建一个空的axios 对象
      // const instance=axios.create({
      //     withCredentials: true,      // 如果发送请求的时候需要带上token 验证之类的也可以写在这个对象里
      //     headers: {
      //         'Content-Type':'multipart/form-data'
      //     }
      // })
      console.log(formData);
      // this.axios.post('http://127.0.0.1:8000/backend/upload/',formData).then(res=>{
      //     if(res.data.code === 200) {
      //         alert(res.data.msg);
      //         this.checkStatus = res.data;
      //     }else if(res.data.code === 2) {
      //         alert(res.data.msg)
      //     }else{
      //         alert(res.data.msg);
      //     }
      // })

      this.axios({
        url: "selina/",
        method: "post",
        data: formData,
        labelPosition: "left",
        headers: {
          "Content-Type": "multipart/form-data",
        },
      }).then((res) => {
        console.log(res);
      });
      alert("The result will be sent to your email :)");
    },
  },
};
</script>

<style scoped>
.el-select {
  width: 90%;
}

.el-input {
  width: 90%;
}

#annotation {
  overflow-x: hidden;
}

.anno {
  /* background-color: rgb(26, 96, 67); */
  border-radius: 30px;
  border: 3px solid rgb(26, 96, 67);
  padding: 10px;
  width: 100%;
  /* display: flex;
    flex-wrap: wrap;

    flex-direction: column; */
  /* align-items: center; */
  font-weight: bold;
}

strong {
  font-size: 30px;
}

.el-col {
  margin-top: 15px;
  margin-bottom: 15px;
  /* margin-left:5px; */
}

#title {
  display: flex;
  flex-wrap: wrap;

  flex-direction: column;
  justify-content: center;
}

.el-button {
  background-color: rgb(26, 96, 67);
}

#help-button{
    /* opacity: 50%; */
    /* width: 2em; */
      /* background-color: rgb(26, 96, 67); */

}

.el-divider {
  color: rgb(26, 96, 67);
  background-color: rgb(26, 96, 67);
}

#Format {
  border-radius: 20px;
  border: 3px solid rgb(26, 96, 67);
  margin: 10px;
  padding: 10px;
}

.custom-icon {
   /* font-size: 2em; */
}

/* ul {
    margin: 20px;
} */

code {
  font-family: Consolas, "courier new";
  color: crimson;
  background-color: #f1f1f1;
  padding: 2px;
  font-size: 105%;
}

.doclink {
  color: rgb(26, 96, 67);
  font-weight: bold;
}

/* CSS link color (red) */
.doclink:hover {
  color: rgb(3, 138, 81);
}

.help {
    height: 30px;
}
/* .file-btn{
    left:0;
    position: relative;
} */

/* .file {
  position: relative;
  display: inline-block;
  background: #d0eeff;
  border: 1px solid #99d3f5;
  border-radius: 4px;
  padding: 4px 12px;
  overflow: hidden;
  color: #1e88c7;
  text-decoration: none;
  text-indent: 0;
  line-height: 20px;
  font-size: 100px;
   width: 60px;
            height: 20px;
}
.file input {
  position: absolute;
  font-size: 100px;
  right: 0;
  top: 0;
  opacity: 0;
}
.file:hover {
  background: #aadffd;
  border-color: #78c3f3;
  color: #004974;
  text-decoration: none;
} */
</style>

