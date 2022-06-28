<template>
  <div id="gene">
    <el-card>
      <el-row :gutter="20" style="margin-top: 50px">
      </el-row>
      <el-row :gutter="20">
        <div style="margin-top: 15px">
          <el-input
            placeholder="Input gene"
            v-model="inputGene"
            class="input-with-select"
            ref="mark"
            @change="isShow"
          >
            <el-select
              v-model="tissuename"
              slot="prepend"
              placeholder="Select tissue"
              @change="isShow"
            >
              <el-option
                v-for="item in options"
                :key="item.value"
                :label="item.label"
                :value="item.value"
              >
              </el-option>
            </el-select>
            <el-button icon="el-icon-search" slot="append"  @click="getHeatmap()"></el-button>
          </el-input>
        </div>
      </el-row>

      <el-row :gutter="20" style="margin-top: 50px">
        <el-col :span="24">
          <div v-if="show_gene_heatmap">
            <el-col :span="12">
              <img
                class="image3"
                :src="`data:image/png;base64,${img_heatmap}`"
              />
            </el-col>
          </div>
        </el-col>
      </el-row>
    </el-card>
  </div>
</template>

<script>
import allts from "../../public/alltissue.json";
import { getHeatmap } from "@/api/getGene.ts";

export default {
  data() {
    return {
      options: allts,
      tissuename: "",
      inputGene: "",
      select: "",
      show_gene_heatmap: false,
    };
  },
  created() {

  },
  
  methods: {
    isShow() {
      this.show_gene_heatmap = false;
    //   console.log("XXX");
    //   console.log(this.show_gene_heatmap);
    },

    async getHeatmap() {
      //   console.log(this.show_gene_violin);
      try{

      
      console.log(this.inputGene);
      this.$refs['mark'].focus()

      const response = await this.axios({
        method: "post",
        url: "get_heatmap/",
        data: {
          gene: this.inputGene,
          tissuename: this.tissuename,
        },
        headers: {
          "Content-Type": "application/json",
        },
      });
      this.img_heatmap = response["data"];
      this.show_gene_heatmap = true;

      console.log(this.show_gene_heatmap);
      }catch(err){
        this.$alert('Please check the input genes, if you are not sure, you can browse our Datasets page to see related genes.', 'Gene', {
          confirmButtonText: 'OK',
          callback: action => {
            this.$message({
              type: 'info',
              message: `action: ${ action }`
            });
              }
        });
      }
    },
  },
};
</script>

<style  >
.el-input {
  width: 600px;
}

.el-select > .el-input {
  width: 300px;
}

.input-with-select .el-input-group__prepend {
  background-color: #fff;
}
</style>
