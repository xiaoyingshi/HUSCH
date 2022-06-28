<template>
  <div id="info">
    <!-- <router-view v-if="isRouterAlive"></router-view> -->
    <el-row>
      <el-col>
        <div class="grid-content bg-purple-dark"></div>
        <el-card class="box-card">
          <el-row :gutter="20">
            <el-select v-model="input" filterable clearable placeholder="Select Tissue">
              <el-option
                v-for="item in TissueFilter"
                :key="item.value"
                :label="item.label"
                :value="item.value"
              >
              </el-option>
            </el-select>
          </el-row>
          <el-table
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
            stripe
            border
            :default-sort="{ prop: 'Tissue', order: 'ascending' }"
            @sort-change="sortChange"
            :row-class-name="tableRowClassName"
            @filter-change="filterChange"
            @row-click="onRowClick"
          >
         

            <template slot="empty">
              <div class="noData">No Results :(</div>
            </template>

            <el-table-column
              prop="Tissue"
              label="Tissue"
              align="center"
              min-width="15%"
              sortable="custom"
              column-key="Tissue"
            >
            </el-table-column>

            <el-table-column
              prop="Dataset"
              label="Dataset"
              align="center"
              min-width="30%"
              @click=""
              sortable="custom"
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
              sortable="custom"
            >
            </el-table-column>
            <el-table-column
              prop="PMID"
              label="PMID"
              align="center"
              min-width="10%"
              sortable="custom"
            >
              <template slot-scope="scope">
                <a class="doclink" :href="scope.row.PMID" target="_blank">{{
                  scope.row.PMIDid
                }}</a>
              </template>
            </el-table-column>
            <el-table-column
              prop="GSE"
              label="Source"
              align="center"
              min-width="15%"
              sortable="custom"
            >
              <template slot-scope="scope">
                <a :href="scope.row.GSE" target="_blank" class="doclink">{{
                  scope.row.source
                }}</a>
              </template>
            </el-table-column>
            <el-table-column
              prop="Stage"
              label="Stage"
              align="center"
              min-width="10%"
              sortable="custom"
            >
            </el-table-column>
            <el-table-column
              prop="Platform"
              label="Platform"
              align="center"
              min-width="15%"
              sortable="custom"
            >
            </el-table-column>
            <el-table-column
              prop="CellNumber"
              label="CellNumber"
              align="center"
              min-width="15%"
              sortable="custom"
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
              :total="
                data_list.filter(
                  (data) => !input || data.Tissue.includes(input)
                ).length
              "
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
// import SelectHeader from "@/components/SelectHeader";

export default {
  name: "info",

  data() {
    return {
      itemList: [],
      data_list: data_json,
      res: "",
      input: "",
      input_tissue: "",
      hover: false,
      DatasetName: "",
      currentPage: 1,
      pageSize: 10,
      isRouterAlive: true,
      TissueFilter: [
                { value: 'Adrenal-Gland', label: 'Adrenal-Gland' },
                { value: 'Airway-Epithelium', label: 'Airway-Epithelium' },
                { value: 'Artery', label: 'Artery' },
                { value: 'Bladder', label: 'Bladder' },
                { value: 'Blood', label: 'Blood' },
                { value: 'Bone-Marrow', label: 'Bone-Marrow' },
                { value: 'Brain', label: 'Brain' },
                { value: 'Breast', label: 'Breast' },
                { value: 'Choroid', label: 'Choroid' },
                { value: 'Decidua', label: 'Decidua' },
                { value: 'Esophagus', label: 'Esophagus' },
                { value: 'Eye', label: 'Eye' },
                { value: 'Fallopian', label: 'Fallopian' },
                { value: 'Gall-Bladder', label: 'Gall-Bladder' },
                { value: 'Heart', label: 'Heart' },
                { value: 'Intestine', label: 'Intestine' },
                { value: 'Kidney', label: 'Kidney' },
                { value: 'Liver', label: 'Liver' },
                { value: 'Lung', label: 'Lung' },
                { value: 'Muscle', label: 'Muscle' },
                { value: 'Nose', label: 'Nose' },
                { value: 'Ovary', label: 'Ovary' },
                { value: 'Pancreas', label: 'Pancreas' },
                { value: 'Peritoneum', label: 'Peritoneum' },
                { value: 'Placenta', label: 'Placenta' },
                { value: 'Pleura', label: 'Pleura' },
                { value: 'Prostate', label: 'Prostate' },
                { value: 'Skin', label: 'Skin' },
                { value: 'Spleen', label: 'Spleen' },
                { value: 'Stomach', label: 'Stomach' },
                { value: 'Testis', label: 'Testis' },
                { value: 'Thyroid', label: 'Thyroid' },
                { value: 'Ureter', label: 'Ureter' },
                { value: 'Uterus', label: 'Uterus' },
                { value: 'Visceral-Adipose', label: 'Visceral-Adipose' },
 



            ]
    };
  },
  change(e) {
    this.$forceUpdate();
  },
  created() {
    //   this.TissueFilter = getFilter();
  },
  computed: {},

  methods: {
    //   renderHeaderOne(h, { column, $index }) {
    //   return (
    //     <span>
    //       {column.label}
    //       <el-popover placement='bottom' width='200' height='200' trigger='click' v-model={this.visible}>
    //         <span slot='reference'>
    //           <i class='el-icon-search' style={this.input ? {'color': 'red'} : {'color': 'blue'}}></i>
    //         </span>
    //           <el-select v-model="value" filterable placeholder="请选择">
    //           <el-option v-for={item in this.options} :key={item.value} :label="item.label" :value="item.value"></el-option>
    // </el-select>
    //         // <el-input size='small' v-model={this.input} placeholder=''></el-input>
    //       </el-popover>
    //     </span>
    //   )
    // },
    renderHeaderOne(h) {
      //下拉框选项
      //   let filters = [
      //     { text: "全部", value: "全部" },
      //     { text: "INFO", value: "INFO" },
      //     { text: "WARN", value: "WARN" },
      //     { text: "ERROR", value: "ERROR" },
      //   ];

      var data = JSON.parse(JSON.stringify(data_json));
      //   console.log(data.length);
      var TissueDict = [];
      for (var i = 0; i < data.length; i++) {
        TissueDict[i] = data[i].Tissue;
      }
      this.TissueList = [...new Set(TissueDict)].sort();
      //   console.log(this.TissueList);
      var TissueFilter = [];
      for (var i in this.TissueList) {
        var tmp_dict = {};
        this.ts = this.TissueList[i];
        tmp_dict["value"] = this.ts;
        tmp_dict["label"] = this.ts;
        TissueFilter.push(tmp_dict);
      }
      //   let filters=TissueFilter

      //下拉框内容包裹在一个div里面
      return h(
        "div",
        // {
        //   style: {
        //     height: "56px",
        //     width: "60px",
        //   },
        // },
        [
          h(
            "span",
            {
              //div里面有一个文字提示：下拉框所属内容
              style: {},
              class: "level-font-class",
            },
            "Tissue"
          ),
          h(
            "el-select",
            {
              //el-select实现下拉框
              on: {
                input: (value) => {
                  //随着下拉框的不同，文字框里的内容在边
                  this.logLevel = value;
                },
              },
              props: {
                value: this.logLevel, //文字框的内容取决于这个value，如果value不存在，会报错
              },
            },
            [
              //下拉框里面填充选项，通过filters遍历map，为每一个选项赋值。
              TissueFilter.map((item) => {
                return h("el-option", {
                  props: {
                    value: item.value,
                    label: item.text,
                  },
                });
              }),
            ]
          ),
        ]
      );
    },

    getFilter() {
      var data = JSON.parse(JSON.stringify(data_json));
      //   console.log(data.length);
      var TissueDict = [];
      for (var i = 0; i < data.length; i++) {
        TissueDict[i] = data[i].Tissue;
      }
      this.TissueList = [...new Set(TissueDict)].sort();
      //   console.log(this.TissueList);
      var TissueFilter = [];
      for (var i in this.TissueList) {
        var tmp_dict = {};
        this.ts = this.TissueList[i];
        tmp_dict["text"] = this.ts;
        tmp_dict["value"] = this.ts;
        TissueFilter.push(tmp_dict);
      }
      //   console.log(TissueFilter);
      return TissueFilter;
    },

    filterHandler(value, row, column) {
      const property = column["property"];
      return row[property] === value;
    },

    fnFilterChangeInit(filter) {},

    filterChange(filterObj) {
      //   console.log(JSON.stringify(filterObj.Tissue[0]));
      this.input = JSON.stringify(filterObj.Tissue[0]).replace(
        /^"(.*)"$/,
        "$1"
      );
    },

    sortChange({ prop, order }) {
      //   console.log("sort dataset");
      this.data_list.sort(this.compare(prop, order));
    },
    compare(propertyName, sort) {
      return function (obj1, obj2) {
        var value1 = obj1[propertyName];
        var value2 = obj2[propertyName];
        if (typeof value1 === "string" && typeof value2 === "string") {
          const res = value1.localeCompare(value2, "zh");
          return sort === "ascending" ? res : -res;
        } else {
          if (value1 <= value2) {
            return sort === "ascending" ? -1 : 1;
          } else if (value1 > value2) {
            return sort === "ascending" ? 1 : -1;
          }
        }
      };
    },

    filData() {},
    // 渲染**的tableheader
    renderSpecNameHeader(createElement, { column, $index }) {
      const self = this;
      // 该列的绑定数据
      // console.log(column.label);
      // 列号
      // console.log($index);
      return createElement(
        "div",
        {
          style: "display:inline-flex;",
        },
        [
          createElement("div", {
            domProps: {
              innerHTML: column.label,
            },
          }),
          createElement(SelectHeader, {
            style: "cursor: pointer;",
            // 组件 prop
            props: {
              type: column.property,
              options: self.specIdOptions, // 下拉框选项
              defaultValue: self.examinerFieldChname, // 默认值
              defaultProps: {
                value: "examinerFieldName",
                label: "examinerFieldChname",
              },
            },

            on: {
              selectChange: self.selectChange,
              resetChange: self.resetChange,
              // click: this.clickHandler
            },
            // 仅用于组件，用于监听原生事件，而不是组件内部使用
            // `vm.$emit` 触发的事件。
            nativeOn: {
              // click: this.nativeClickHandler
            },
          }),
        ]
      );
    },
    // 选择框回调
    selectChange(data) {
      console.log("回调", data);
      // 自定义筛选框返回数据进行过滤添加到tableData数组中
      const type = data["type"];
      const value = data["value"];
      this.rules[type] = value;
      if (value !== "" && type !== "") {
        this.tableDataDeal = this.tableDataDeal.filter((item) =>
          item[type]
            .toString()
            .toLowerCase()
            .includes(value.toString().toLowerCase())
        );
      }
    },
    // 重置回调
    resetChange(data) {
      console.log("重置回调", data);
      delete this.rules[data["type"]];
      var tmpData = this.tableDataDealCopy;
      for (const key in this.rules) {
        tmpData = tmpData.filter((item) =>
          item[key]
            .toString()
            .toLowerCase()
            .includes(this.rules[key].toString().toLowerCase())
        );
      }
      this.tableDataDeal = tmpData;
    },

    indexNum(index) {
      return index + (currentPage - 1) * 10 + 1;
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
      //   console.log(this.DatasetName);
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

      console.log(this.$router);
      // this.$router.go(0);
    },
  },
};
</script>


<style lang="scss" scoped >
@import "@/styles/info.scss";
</style>