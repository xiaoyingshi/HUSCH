<template>
  <div id="info">
    <el-row>
      <el-card>
        
        <div id="tiss" type="primary" round>
          <img
                      height="70px"
                      class="image"
                      v-bind:src="
                        require('@/assets/' +
                          input +
                          '.png')
                      "
                    />
                    <p>{{input}}</p>
          </div>
        
      
      <div>
        <h1>Celltype</h1>
        <!-- <h3>Celltype-Dataset</h3>
        <img  width="95%" src="@/assets/sum_ct_eg.png" alt=""> -->
         <h3>Tissue-level</h3>
        <!-- <img  width="95%" src="@/assets/skin_ct.png" alt=""> -->
        <div id='clu'>
                    <div style="width:48%">
                      <img
                      
                      class="image"
                      src="@/assets/skin_ct.png"
                      style="max-width: 100%; max-height: 100%"
                    />
                    </div>
                    <div style="width:48%">
                      <img
                      
                      class="image"
                      src="@/assets/skin_ct.png"
                      style="max-width: 100%; max-height: 100%"
                    />
                    </div>
                    
                    <div id="harmony">
                      <div>
                        <el-select v-model="ct_value" placeholder="Select celltype" @change="getPltUrl">
                    <el-option
                      v-for="item in pic_type"
                      :key="item.value"
                      :label="item.label"
                      :value="item.value">
                    </el-option>
                  </el-select>
                      </div>
                      
                    <img
                    width="48%"
                      class="image"
                      v-if="plt_show"
                      src="@/assets/skin_ct.png"
                    
                    />
                    </div>
                    
                  
                  </div>
      </div>
      <div>
        <h1>Marker</h1>
        <el-row>
          <el-col>
            <el-table
        :data="tableData"
        stripe
        style="width: 100%"
        height='300'
        >
        <el-table-column
        prop='Celltype'
        label='Celltype'
        width="180"
        >

        </el-table-column>
        <el-table-column
        prop='Marker'
        label='Marker'
        >

        </el-table-column>

        </el-table>

          </el-col>
        </el-row>
        
      </div>
      </el-card>
      
    </el-row>
    <!-- <router-view v-if="isRouterAlive"></router-view> -->
    <el-row>
      <!-- <el-col :span="15"> -->
      <el-col>
        <div class="grid-content bg-purple-dark"></div>
        <el-card class="box-card">
          <h3>Dataset</h3>
          <el-table
            v-model="input"
            :data="
              data_list
                .filter(
                  (data) =>
                    !input ||
                    data.Tissue.toLowerCase().includes(input.toLowerCase()) ||
                    data.Dataset.toLowerCase().includes(input.toLowerCase())
                )
                .slice((currentPage - 1) * pageSize, currentPage * pageSize)
            "
            border
            :default-sort="{ prop: 'Tissue', order: 'ascending' }"
            @sort-change
            :row-class-name="tableRowClassName"
            @row-click="onRowClick"
          >
            <template slot="empty">
              <div class="noData">No Results :(</div>
            </template>

            <el-table-column type="index" label="#"></el-table-column>
            <el-table-column
              prop="Tissue"
              label="Tissue"
              align="center"
              min-width="15%"
              sortable
            >
            </el-table-column>
            <el-table-column
              prop="Dataset"
              label="Dataset"
              align="center"
              min-width="30%"
              @click=""
              sortable
            >
              <template slot-scope="scope">
                <el-popover placement="right" width="200" trigger="hover">
                  <span>
                    <img
                      :src="
                        require(`@/assets/pictures/${scope.row.Dataset}_cluster.png`)
                      "
                      style="max-height: 500px; max-width: 500px"
                    />
                    <img
                      :src="
                        require(`@/assets/pictures/${scope.row.Dataset}_assign.png`)
                      "
                      style="max-height: 500px; max-width: 500px"
                    />
                  </span>
                  <div slot="reference">
                    {{ scope.row.Dataset }}
                  </div>
                </el-popover>
              </template>
            </el-table-column>
            <el-table-column
              prop="Year"
              label="Year"
              align="center"
              min-width="10%"
              sortable
            >
            </el-table-column>
            <el-table-column
              prop="PMID"
              label="PMID"
              align="center"
              min-width="10%"
              sortable
            >
              <template slot-scope="scope">
                <a
                  :href="'https://pubmed.ncbi.nlm.nih.gov/' + scope.row.PMID"
                  target="_blank"
                  class="buttonText"
                  >{{ scope.row.PMID }}</a
                >
              </template>
            </el-table-column>
            <el-table-column
              prop="GSE"
              label="GSE"
              align="center"
              min-width="15%"
              sortable
            >
            </el-table-column>
            <el-table-column
              prop="Stage"
              label="Stage"
              align="center"
              min-width="10%"
              sortable
            >
            </el-table-column>
            <el-table-column
              prop="Platform"
              label="Platform"
              align="center"
              min-width="15%"
              sortable
            >
            </el-table-column>
          </el-table>

          <div class="block" style="text-align: center">
            <el-pagination
              @size-change="handleSizeChange"
              @current-change="handleCurrentChange"
              :current-page.sync="currentPage"
              :page-size="10"
              layout="prev, pager, next, jumper"
              :total="data_list.length"
            >
            </el-pagination>
          </div>
        </el-card>
      </el-col>
    </el-row>
  </div>
</template>

<script>
import axios from "axios";
import data_json from "@/assets/data.json"; // 引用

export default {
  name: "info",
  
  
  data() {
    return {
      itemList: [],
      data_list: data_json,
      res: "",
      input:this.$route.params.tissue,
      input_tissue: "",
      hover: false,
      DatasetName: "",
      currentPage: 1,
      pageSize: 10,
      isRouterAlive: true,
      level1:"",
      level3:"",
      pic_type:[
        {value:'Immune',label:'Immune'},
        {value:'Stromal',label:'Stromal'},
        {value:'Tissue-specific',label:'Tissue-specific'},
        {value:'Endothelial',label:'Endothelial'}
      ],
      ct_value:"",
      plt_show:false,
      tableData: [
        { Celltype : "ORS B" , Marker : "GJB6;KRT16;KRT6A;KRT6B"},
{ Celltype : "ORS SB" , Marker : "GJB6"},
{ Celltype : "IRS H/H" , Marker : "FABP9"},
{ Celltype : "IFE Spinous" , Marker : "SPINK5;CALML5;KRTDAP;KRT10;KRT1"},
{ Celltype : "IFE Basal" , Marker : "COMP;MGP;CXCL14;HAA1;CHISL1"},
{ Celltype : "ORS CL" , Marker : "KRT75;TM4SF1"},
{ Celltype : "IFE Granular" , Marker : "SPINK5;CALML5;CST6;KRT16;KRT6A;KRT6B"},
{ Celltype : "Matrix/cortex/medulla" , Marker : "KRT85;KRT35"},
{ Celltype : "IFE Mitotic" , Marker : "HMGB2;HIST1H4C"},
{ Celltype : "Isthmus" , Marker : "CYR61;EPCAM"},
{ Celltype : "IRS/cuticle" , Marker : "KRT28;KRT73;FABP9;TCHH"},
{ Celltype : "Mesenchymal" , Marker : "CYR61;DCN;CFD;MGP"},
{ Celltype : "Melanocyte" , Marker : "MLANA;PMEL;PMEL;MLANA;TYRP1;DCT"},
{ Celltype : "Infundibulum" , Marker : "RCAN1;KRT15"},
{ Celltype : "Sebaceous/apocrine" , Marker : "DCD;MUCL1"},
{ Celltype : "Bulge" , Marker : "RCAN1;CXCL14;KRT15"},
{ Celltype : "Lower Bulge" , Marker : "COMP;MGP;KRT15"},
{ Celltype : "Endothelial" , Marker : "TM4SF1;SPARCL1;IFI27"},
{ Celltype : "Langerhans" , Marker : "SRGN;CD74;HLA-DRA"},
{ Celltype : "Immune/Ts" , Marker : "SRGN;CD69"},
{ Celltype : "CD4T" , Marker : "CD4;CD3;CD4;CD3"},
{ Celltype : "CD8T" , Marker : "CD8A;CD8B;CD3;CD8A;CD8B;CD3"},
{ Celltype : "Pericyte" , Marker : "ACTA2;RGS5;PDGFRB"},
{ Celltype : "Keratinocyte" , Marker : "KRT1;KRT10;SBSN;KRTDAP"},
{ Celltype : "Fibroblast" , Marker : "PDGFRA;LUM;DCN;VIM;COL1A2"},
{ Celltype : "Macrophages/DC" , Marker : "AIF1;LYZ;HLA-DRA;CD68;ITGAX"},
{ Celltype : "T" , Marker : "CD3D;CD3G;CD3E;LCK"},
{ Celltype : "Vascular Endothelial" , Marker : "SELE;CLDN5;VWF;CDH5"},
{ Celltype : "Lymphatic Endothelial" , Marker : "PROX1;CLDN5;LYVE1"},
{ Celltype : "Erythrocyte" , Marker : "HBA1;HBA2;HBB"},
{ Celltype : "Epidermal Stem" , Marker : "KRT5;KRT14;TP63;ITGB1;ITGA6"}
      ]
    };
  },
  change(e) {
    this.$forceUpdate();
  },
  created() {
    // this.getData();
    // this.getDataList();
    function isFileExisted(file) {
    return new Promise((resolve, reject) => {
      fs.access(file, (err) => {
        if (err) {
            this.level1='@/assets/integrate/Adipose_harmony_level1.png'
            this.level3='@/assets/integrate/Adipose_harmony_level3.png'
          resolve(false);//"不存在"
        } else {
            console.log('cunzai')
          resolve(true);//"存在"
        }
      })
    })
};
this.le


    
  },
  computed: {},
  
  methods: {
    // async getDataList() {
    //   const { data: res } = await this.axios.get("./data.json", {
    //     params: this.queryInfo,
    //   });
    //   // console.log(res);
    //   this.data_list = res;
    //   // this.total = res.data.total;
    //   // console.log(this.total);
    // },
    getPltUrl(){
      this.plt_url=this.DatasetName+'_harmony_'+this.ct_value
      this.plt_show=true
    },
    tableRowClassName({ row, rowIndex }) {
      //Put the index of each row into row
      row.index = rowIndex;
    },
    handleSizeChange(val) {
      this.currentPage = 1;
      this.currentPage = val;
      this.pageSize = val;
      // console.log(`每页 ${val} 条`);
    },
    handleCurrentChange(val) {
      this.currentPage = val;
    },

    onRowClick(row) {
      //Click on the row to get the index
      this.DatasetName = row.Dataset;
      console.log(this.DatasetName);
      // need return use push
      // this.$router.push('/detail/');
      this.$router.push(
        "/detail/" + this.DatasetName
        // {
        // name: "detail",
        // params: {
        //   DatasetName: this.DatasetName,
        // },

        // }
        
      );
      
      
      console.log('hahah');
      console.log(this.$router);
      // this.$router.go(0);
    },
  },
};
</script>


<style lang="scss" >
@import "@/styles/info.scss";
#tiss{
  font-size: 60px;
  display: flex;
  flex-direction: row;
  align-items: center;
  
}
#clu{
  display: flex;
  flex-direction: row;
}
#harmony{
  display: flex;
  flex-direction: column;
}
</style>