<template>

  <div id="info" >
    
    <el-row>
      <el-col>
        <div class="grid-content bg-purple-dark"></div>
        <el-tabs type="border-card">
          <el-tab-pane
            ><span slot="label" class="tab-header">Overview</span>
            <el-row>
              <el-col :span="24">
                <el-card>
                  <div v-cloak>
                    <h2 style="text-align: center">{{ DatasetName }}</h2>
                  </div>
                  <div id='clu'>
                    <el-col :span="12" >
                    <img
                      class="image"
                      v-bind:src="
                        require('@/assets/pictures/' +
                          DatasetName +
                          '_cluster.png')
                      "
                    />
                  </el-col>
                  <el-col :span="12">
                    <img
                      class="image"
                      v-bind:src="
                        require('@/assets/pictures/' +
                          DatasetName +
                          '_assign.png')
                      "
                    />
                  </el-col>
                  </div>
                  
                </el-card>
              </el-col>
            </el-row>
            <el-row>
              <el-col :span="24">
                <el-card>
                
                <h4 style="text-align: center">
                  Cell-type statistics
                </h4>
                <div id='clu'>
                    <el-col :span="12" >
                    <img
                      class="image"
                      v-bind:src="
                        require('@/assets/pictures/' +
                          DatasetName +
                          '_Pie.png')
                      "
                    />
                  </el-col>
                  <el-col :span="12">
                    <img
                      class="image"
                      v-bind:src="
                        require('@/assets/pictures/' +
                          DatasetName +
                          '_ProBar.png')
                      "
                    />
                  </el-col>
                  </div>
                
              </el-card>
              </el-col>
              
            </el-row>
            <el-row>
              <el-card shadow="always">
                <h4 style="text-align: center">
                  Top differential genes for each cell type
                </h4>
                <el-table
                  ref="filterTable"
                  :data="json_list"
                  stripe
                  style="width: 100%"
                  max-height="250"
                  v-loading="t_loading"
                  
                >
                  <el-table-column
                    
                    prop="cluster"
                    label="Cell Type"
                    
                    align="center"
                    min-width="10%"
                    column-key="cluster"
                    :filters="ct_list"
                    :filter-method="filterHandler"
                    
                  >
                  </el-table-column>

                  <el-table-column
                    prop="gene"
                    label="Gene"
                    align="center"
                    min-width="15%"
                  >
                  </el-table-column>
                  <el-table-column
                    prop="p_val_adj"
                    label="Adjusted p-value"
                    align="center"
                    sortable
                    min-width="15%"
                  >
                  </el-table-column>
                  <el-table-column
                    prop="avg_logFC"
                    label="log2FC"
                    align="center"
                    sortable
                    min-width="30%"
                  >
                  </el-table-column>
                </el-table>

                <!-- <div class="block" style="text-align: center">
                  <el-pagination
                    @size-change="handleSizeChange"
                    @current-change="handleCurrentChange"
                    :current-page.sync="currentPage"
                    :page-size="10"
                    layout="prev, pager, next, jumper"
                    :total="json_list.length"
                  >
                  </el-pagination>
                </div> -->
              </el-card>
            </el-row>
          </el-tab-pane>
          <el-tab-pane
            ><span slot="label" class="tab-header">Gene</span>
            <el-card>
              <el-row>
                <el-col :span="6"> Individual gene exploration </el-col>
              </el-row>
              <el-row>
                <!-- <el-input
                  class="inline-input"
                  v-model.lazy="inputGene"
                  :fetch-suggestions="querySearch"
                  placeholder="Input gene"
                  @select="handleSelect"
                  style="width:25%"
                > -->
                <el-select
    v-model="inputGene"
    multiple
    filterable
    remote
    reserve-keyword
    placeholder="Please input Gene"
    :remote-method="remoteMethod"
    :loading="loading">
                  <el-option
                    v-for="item in optionss"
                    :key="item.value"
                    :label="item.label"
                    :value="item.value">
                  </el-option>
                </el-select>
                  <el-button
                    
                    icon="el-icon-search"
                    @click="searchGene()"
                  >Search</el-button>
                
              </el-row>
              <el-row :gutter="20" style="margin-top: 50px">
                <el-col :span="24">
                  <div v-cloak>
                    <h2 style="text-align: center">{{ DatasetName }}</h2>
                  </div>
                  <el-col :span="12">
                    <img
                      class="image2"
                      v-bind:src="
                        require('@/assets/pictures/' +
                          DatasetName +
                          '_assign.png')
                      "
                      width="100px"
                      height="300px"
                    />
                  </el-col>
                  <div v-if="show_gene_umap">
                    <el-col :span="12">
                      <img
                        class="image3"
                        :src="`data:image/png;base64,${img}`"
                        
                      />
                    </el-col>
                  </div>
                </el-col>
              </el-row>
            </el-card>

            <el-row>
              <el-card
                ><el-row style="margin-bottom: 50px">
                  <el-row>
                    <el-col :span="6"> Comparison </el-col>
                  </el-row>
                  <!-- <el-col :span="6"> Comparison </el-col> -->
                  <el-select
                    v-model="compareType"
                    placeholder="Select"
                    @change="isShow()"
                    value-key="value"
                  >
                    <el-option
                      v-for="item in options"
                      :key="item.value"
                      :label="item.label"
                      :value="item.value"
                    >
                    </el-option>
                  </el-select>
                  <el-button type="info" plain @click="getViolin()"
                    >Update</el-button
                  >
                </el-row>
                <el-row>
                  <div>
                    <el-col :span="24">
                      <img
                        v-if="show_gene_violin"
                        class="image4"
                        :src="`data:image/png;base64,${img_viloin}`"
                        width="100%"
                        
                      />
                    </el-col>
                  </div>
                </el-row>
              </el-card>
            </el-row>
          </el-tab-pane>

          <el-tab-pane>
            <span slot="label" class="tab-header">GSEA</span>
            <el-card>
              <el-tabs type="border-card">
                <el-tab-pane>
                  <span slot="label" class="tab-header">KEGG</span>
                  
                    <el-image 
                      
                      v-bind:src="
                        require('@/assets/pictures/' +
                          DatasetName +
                          '_kegg_UP_heatmap.png')
                      "
                      :preview-src-list="[require('@/assets/pictures/'+DatasetName +'_kegg_UP_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_kegg_DOWN_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_hallmark_UP_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_hallmark_DOWN_heatmap.png')]"
                    >
                    </el-image>  
                 
                    <el-image
                      class="image"
                      v-bind:src="
                        require('@/assets/pictures/' +
                          DatasetName +
                          '_kegg_DOWN_heatmap.png')
                      "
                      :preview-src-list="[require('@/assets/pictures/'+DatasetName +'_kegg_DOWN_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_kegg_UP_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_hallmark_UP_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_hallmark_DOWN_heatmap.png')]"
                    >
                  </el-image>
                  
                </el-tab-pane>
                <el-tab-pane>
                  <span slot="label" class="tab-header">Hallmark</span>
                  
                    
                    <el-image
                      class="image"
                      v-bind:src="
                        require('@/assets/pictures/' +
                          DatasetName +
                          '_hallmark_UP_heatmap.png')
                      "
                      :preview-src-list="[require('@/assets/pictures/'+DatasetName +'_hallmark_UP_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_hallmark_DOWN_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_kegg_UP_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_kegg_DOWN_heatmap.png')]"
                    >
                  </el-image>
                    <el-image
                      class="image"
                      v-bind:src="
                        require('@/assets/pictures/' +
                          DatasetName +
                          '_hallmark_DOWN_heatmap.png')
                      "
                      :preview-src-list="[require('@/assets/pictures/'+DatasetName +'_hallmark_DOWN_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_hallmark_UP_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_kegg_UP_heatmap.png'),require('@/assets/pictures/'+DatasetName +'_kegg_DOWN_heatmap.png')]"
                    >
                  </el-image>
                </el-tab-pane>
                
              </el-tabs>
            </el-card>
          </el-tab-pane>
          <el-tab-pane>
            <span slot="label" class="tab-header">Cell-Cell Interaction</span>
            
              <el-row :span="24">
                <el-card>
                   <!-- <el-col :span="12"> -->
                     <div id='cc'>
                       <el-col :span="12">
                        <img
                      class="image"
                      
                      v-bind:src="
                        require('@/assets/pictures/' +
                           DatasetName+
                          '_log_CCI.png')
                      "
                    />
                    </el-col>
                    <el-col :span="12">
                    <img
                      class="image"
                      
                      v-bind:src="
                        require('@/assets/pictures/' +
                           DatasetName+
                          '_circle_CCI.png')
                      "
                    />
                    </el-col>
                    </div>
                </el-card>
              <el-select v-model="value" placeholder="Celltype" @change="getPltUrl">
              <el-option
                v-for="item in CT"
                :key="item"
                :label="item"
                :value="item">
              </el-option>
              
          </el-select>
          <div v-if="warning" id="warning">
            <h2>Sorry, this dataset has not been updated yet, the example results are shown as follows : </h2>
            <el-col :span="12">
                        <img
                      class="image"
                      src="@/assets/pictures/CCI_CT_plt/HU_0229_Pancreas_GSE85241_Beta_in.png"
                      
                    />
                    </el-col>
                    <el-col :span="12">
                    <img
                      class="image"
                      src="@/assets/pictures/CCI_CT_plt/HU_0229_Pancreas_GSE85241_Beta_out.png"
                    />
                    </el-col>
          </div>
          <div id='cc' v-if="plt_show">
                       <el-col :span="12">
                        <img
                      class="image"
                      
                      v-bind:src="
                        require('@/assets/pictures/CCI_CT_plt/' +
                           plt_url+
                          '_in.png')
                      "
                    />
                    </el-col>
                    <el-col :span="12">
                    <img
                      class="image"
                      
                      v-bind:src="
                        require('@/assets/pictures/CCI_CT_plt/' +
                           plt_url+
                          '_out.png')
                      "
                    />
                    </el-col>
              </div>
              </el-row>
            
          </el-tab-pane>
          <el-tab-pane>
            <span slot="label" class="tab-header">TF</span>
            <el-card>
               <!-- <el-select v-model="value" placeholder="Celltype" @change="getPltUrl">
              <el-option
                v-for="item in CT"
                :key="item"
                :label="item"
                :value="item">
              </el-option>
              
          </el-select> -->
          <img
                      class="image"
                     
                      v-bind:src="
                        require('../assets/pictures/' +
                          DatasetName+
                          '_TF_Heatmap.png')
                      "
                    />
          
            </el-card>
           
          
            
            
          </el-tab-pane>
          <el-tab-pane>
            <span slot="label" class="tab-header">Download</span>
            
              <el-row :span="24">
                <el-card id="download">
                                <!-- <el-col :span="4" > -->
                                <!-- <a :href="'http://39.101.160.21:8000/backend/download/'+DatasetName+'_EM.txt'"><el-button type="primary" round >Expression Matrix</el-button></a> -->
                                <el-button type="primary" round @click="alert">Expression matrix</el-button>
                                <!-- </el-col :span="4"> -->
                                <!-- <el-col> -->

                                <!-- </el-col> -->
                                <!-- <el-col :span="4" > -->
                                <!-- <a :href="'http://39.101.160.21:8000/backend/download/'+DatasetName+'_Meta.txt'"><el-button type="primary" round @click="download()">Mate Info</el-button></a> -->
                                <el-button type="primary" round @click="alert">Meta file</el-button>
                                <!-- </el-col> -->
                            </el-card>
              </el-row>
          </el-tab-pane>
        </el-tabs>

      </el-col>
    </el-row>
  </div>
</template>


<script lang="ts">
import { searchGene } from "@/api/getGene.ts";
import { getGene } from "@/api/getGene.ts";
import { getViolin } from "@/api/getGene.ts";
import { BIconLayoutTextWindowReverse } from "bootstrap-vue";
import TF_CT  from "@/assets/TF_CT.json"; // 引用
// import data_json from "@/assets/JsonDir/"+this.DatasetName+".json"// 引用

export default {
  props: ["img"],
  
  data() {
    return {
      DatasetName: this.$route.params.DatasetName,
      json_list: [],
      ct_list: [],
      ct_tmp:[],
      data_list: [],
      currentPage: 1,
      pageSize: 10,
      inputGene: [],
      available_genes: "",
      state2: "",
      links: [],
      show_gene_umap: false,
      show_gene_violin: false,
      value:"",
      plt_url:"",
      plt_show:false,
      gene:"",
      optionss: [],
      list: [],
      loading: false,
      t_loading: true,
      warning:false,
      // srcList:[
      //   require('@/assets/pictures/' +this.DatasetName +'_hallmark_UP_heatmap.png'),
      //   // require('@/assets/pictures/' +this.DatasetName +'_hallmark_DOWN_heatmap.png')
      //   'https://fuss10.elemecdn.com/8/27/f01c15bb73e1ef3793e64e6b7bbccjpeg.jpeg',
      //   'https://fuss10.elemecdn.com/1/8e/aeffeb4de74e2fde4bd74fc7b4486jpeg.jpeg'
      // ],
      // srcList1:require('@/assets/pictures/' +this.DatasetName +'_hallmark_UP_heatmap.png'),

      
      

      umap_path: "",
      options: [
        {
          value: "CellType",
          label: "CellType",
        },
        {
          value: "Cluster",
          label: "Cluster",
        },
      ],
      compareType: "",
    };
  },
  watch:{
      $route(){
        //跳转到该页面后需要进行的操作
      }
    },
    // beforeCreate() {
    //   this.DatasetName = this.$route.params.DatasetName;
    //   console.log(this.$route.params.DatasetName);
    //   this.TF_CT=TF_CT;
    //   this.CT=this.TF_CT[this.DatasetName];
    //   console.log(this.CT);    
    //   this.getJsonFile(this.DatasetName);
    // },
  created() {
    this.getDataList();
    
    this.DatasetName = this.$route.params.DatasetName;
    console.log(this.$route.params.DatasetName);
    this.TF_CT=TF_CT;
    this.CT=this.TF_CT[this.DatasetName];
    console.log(this.CT);    
    this.getJsonFile(this.DatasetName);
    if(typeof(this.CT)=='undefined'){
      this.warning=true
    }
    console.log('length')
    // console.log(class(this.CT))
    // console.log(this.CT.length)
    // this.json_list=this.data_list[0]
    // console.log('hi')
    // console.log(this.data_list[1])
    // this.ct_list=this.data_list[1]
    // this.url_EM="http://127.0.0.1:8000/backend/download/"+DatasetName+'_EM.txt',
    // this.url_Mate="http://127.0.0.1:8000/backend/download/"+DatasetName+'_mate.txt',
  },

  methods: {
    alert() {
            this.$alert(
                "During the test phase, the data is only available for download from the internal network.",
                "Warning",
                {
                    confirmButtonText: "OK",
                    callback: (action) => {
                        this.$message({
                            type: "info",
                            message: `action: ${action}`,
                        });
                    },
                    }
            );
        },
    changeActive($event) {
        $event.currentTarget.className = 'animate__animated animate__headShake';
      },
      removeActive($event) {
        $event.currentTarget.className = '';
      },
    
      handleChange(val) {
        console.log(val);
      },
    test(){
      console.log(this.inputGene)
      console.log(this.inputGene[1])
    },
      remoteMethod(query) {
        if (query !== '') {
          this.loading = true;
          setTimeout(() => {
            this.loading = false;
            this.optionss = this.list.filter(item => {
              
              return item.label.toLowerCase()
                .indexOf(query.toLowerCase()) > -1;
            });
          }, 200);
        } else {
          this.optionss = [];
        }},
      filterHandler(value, row, column) {
        const property = column['property'];
        console.log(property)
        console.log(value)
        return row[property] == value;
      },
    isShow() {
      this.show_gene_violin = !this.show_gene_violin;
    },
    async getJsonFile(dataset) {1
      const response = await this.axios({
        method: "post",
        url: "get_json/",
        data: {
          dataset: this.DatasetName,
        },
        headers: {
          "Content-Type": "application/json",
        },
      });
      //   console.log(response["data"]);
      // this.data_list = response["data"]
      this.json_list = response["data"][0];
      this.ct_tmp = response["data"][1];
      for (var ct in this.ct_tmp){
        var tmp={text:this.ct_tmp[ct],value:this.ct_tmp[ct]}
        // console.log(tmp)
        // console.log(tmp.value)
        this.ct_list.push(tmp)
        
      }
      this.t_loading = false
      //   console.log(this.json_list);
      // console.log(this.json_list);
      // console.log('haha')
      // console.log(this.data_list[1])
      // return this.data_list;
      // console.log(this.ct_list)
      // return this.json_list,this.ct_list
    },

    async getDataList() {
      const { data: res } = await this.axios.get(
        "/" + this.DatasetName + ".json"
      );
      this.json_list = res;
      // this.t_loading = false
      // console.log(this.json_list);
    },
    getPltUrl(){
      this.plt_url=this.DatasetName+'_'+this.value
      this.plt_show=true
    },
    // async download(){
    //   var file_name= this.DatasetName+'.txt';
    //   var file_address='/data/dongqing/yuzhiguang/mysite/download/';
    // //   this.axios.get('http://127.0.0.1:8000/backend/download?file_name='+file_name).then(function(response){
    // //     console.log('success')
    // //     console.log(response.data)
    // //     }).catch(function(err){console.log(err)});
    // //   windows.open('http://127.0.0.1:8000/backend/download?file_name='+file_name)
    // let formData = new FormData();
    // formData.append('name',file_name);
    // formData.append('add',file_address);
    // let config = {
    //             headers: {
    //                 'Content-Type':'multipart/form-data'
    //             }
    //         };
    // this.axios.post('http://127.0.0.1:8000/backend/download/',formData).then(res=>{
    //             if(res.data.code === 200) {
    //                 alert(res.data.msg);
    //                 this.checkStatus = res.data;
    //             }else if(res.data.code === 2) {
    //                 alert(res.data.msg)
    //             }else{
    //                 alert(res.data.msg);
    //             }
    //         })
    // },
    
    tableRowClassName({ row, rowIndex }) {
      //Put the index of each row into row
      row.index = rowIndex;
    },
    // handleSizeChange(val) {
    //   this.currentPage = 1;
    //   this.currentPage = val;
    //   this.pageSize = val;
    // },
    // handleCurrentChange(val) {
    //   this.currentPage = val;
    // },
    handleSelect(item) {
      console.log(item);
    },
    // return matched strings
    querySearch(queryString, cb) {
      var available_genes = this.available_genes;
      var results = queryString
        ? available_genes.filter(this.createFilter(queryString))
        : available_genes;
      // 调用 callback 返回建议列表的数据
      console.log(results);
      //   console.log(typeof available_genes);
      cb(results);
    },

    getImage(search) {
      console.log(this.umap_path);
      return require(`@/assets/media/${this.umap_path}`);
    },

    createFilter(queryString) {
      return (available_gene) => {
        return (
          available_gene.gene
            .toLowerCase()
            .indexOf(queryString.toLowerCase()) === 0
        );
      };
    },

    async getViolin(gene, dataset, compare) {
      console.log(this.show_gene_violin);
      this.inputGene = this.inputGene.filter(item => {
              
              return item.toUpperCase()
              });
      console.log(this.inputGene);

      const response = await this.axios({
        method: "post",
        url: "violin_gene/",
        data: {
          gene: this.inputGene,
          dataset: this.DatasetName,
          compare: this.compareType,
        },
        headers: {
          "Content-Type": "application/json",
        },
      });
      this.img_viloin = response["data"];
      this.show_gene_violin = true;

      console.log(this.compareType);
    },

    async searchGene(gene, dataset) {
      this.show_gene_umap = true;
      // this.inputGene = this.inputGene.toUpperCase()
      this.inputGene=this.inputGene.filter(item => {
              
              return item.toUpperCase()
              });

      console.log(this.show_gene_umap);
      console.log(this.inputGene);
      const response = await this.axios({
        method: "post",
        url: "detail_gene/",
        data: {
          gene: this.inputGene,
          dataset: this.DatasetName,
        },
        headers: {
          "Content-Type": "application/json",
        },
      });
      this.img = response["data"];
    },

    async getGene(dataset) {
      const response = await this.axios({
        method: "post",
        url: "get_gene/",
        data: {
          dataset: this.DatasetName,
        },
        headers: {
          "Content-Type": "application/json",
        },
      });
        // console.log(response["data"]);

      this.available_genes = response["data"];
      
      this.list = this.available_genes.map(item => {
        return { value: `${item}`, label: `${item}` };
      });
      console.log(this.list)
      return this.available_genes;
    },
  },

  mounted() {
    // this.available_genes = this.loadAllgene();
    this.available_genes = this.getGene();
    
    
  },
};
</script>


<style lang="scss" scoped>
@import "../styles/detail.scss";
#clu{
  display: flex;
  flex-direction: row;
}
#cc{
  display: flex;
  flex-direction: row;
}
#download{
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: center;
}
#pic_gsea{
  display: flex;
  flex-direction: row;
  #pie{
    width: 50%;
  }
}
template{
  font-size: 40px;
}
.el-collapse-item__header{
  font-size: 30px !important;
}
.el-collapse-item >>> .el-collapse-item__header{
  font-size: 30px !important;
  color:rgb(39, 196, 145) ;
}
.el-card ::v-deep .el-card_header{
  font-size: 40px;
  background-color:rgb(39, 196, 145) ;
  color: white;
}
.analysis_head{
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items:center;
  color:rgb(39, 196, 145) ;
  // width: 60%;
  margin-bottom: 40px;
}
#warning{
  margin: 30px;
}



</style>