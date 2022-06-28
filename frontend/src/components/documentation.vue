<template>
    <div id="documentation">

        <el-tabs type="border-card">

            <el-tab-pane label="Documentation">
                <div>
                    <div id='doc'>
                        <h1 id="husch">HUSCH</h1>
                        <h2 id="introduction">Introduction</h2><br>
                        <p>Human Universal Single Cell Hub (HUSCH) is a scRNA-seq database focusing on human universal
                            cell at single-cell resolution, consisting of 136 datasets from 35 human tissues covering
                            7 different platforms. HUSCH provides detailed cell-type annotation, online analysis
                            function and online cell-type annotation function based on <a class="doclink"
                                href="https://github.com/SELINA-team/SELINA.py" target="_blank">SELINA</a>.</p><br>
                        <h2 id="data-collection-and-processing">Data Collection and Processing</h2><br>
                        <p>We collected human normal tissue scRNA-seq studies automatically crawled from public
                            databases. For each collected dataset, a uniform analysis pipeline was adopted to perform
                            quality control and clustering. We used gene markers the original studies provided to
                            automatically annotated them based on the strategy used in <a class="doclink"
                                href="https://github.com/liulab-dfci/MAESTRO" target="_blank">MAESTRO</a>. After
                            processing, we
                            curated the cell-type annotation of all datasets at two levels: major-lineage and
                            minor-lineage. The curation makes the gene expression in different cell types comparable
                            across all datasets within each tissue.</p><br>
                        <h2 id="function-of-husch">Functions of HUSCH</h2><br>
                        <p>Based on the unified data processing, HUSCH presents the analysis results in a user-friendly
                            interface for public accessing, which allows researchers to gain a quick insight into the
                            expression of genes of interest at the single-cell level.</p>
                        <p>HUSCH provides online cell-type annotation function based on <a class="doclink"
                                href="https://github.com/SELINA-team/SELINA.py" target="_blank">SELINA</a>, a deep
                            learning-based
                            framework for single cell assignment with multiple references.</p><br>
                        <h2 id="how-to-start">How to Start</h2><br>
                        <p>In the home page, you can select organ they are interested in, then a summary of the tissue
                            will be shown, which includes dataset basic information such as unique dataset name, paper
                            published year, paper PMID, sequencing platforms and so on. You can also explore the
                            specific dataset in
                            Datasets page.</p>
                        <img class="selinalogo" src="../assets/tmp_home.png"><br>

                        <h2 id="specific-dataset">Explore Specific Dataset</h2><br>
                        <p>There are 6 modules: Overview, Gene, GSEA, CCI, TF, Download.</p><br>
                        <p><b>Overview:</b> The first two UMAP show the Cluster and Annotation of this dataset. Then
                            the pie
                            chart shows the proportion of each cell type. HUSCH also list the top differential genes for
                            each cell type in a table.</p><br>
                        <p><b>Gene:</b> This module includes Individual gene exploration and Comparison.</p><br>
                        <img class="selinalogo" src="../assets/doc_gene.png">

                        <p><b>Individual gene exploration:</b> You can enter your gene of interest into the search
                            bar, the
                            corresponding results will be presented on the right by the feature plot.</p><br>
                        <p>If you want to get more detailed information, please select Cell type or Cluster in the
                            checkbox in the Comparison below. The expression of the target gene in cell-type or cluster
                            will be presented by violin plot.</p><br>
                        <p><b>GSEA:</b> In this module, we provide the gene set enrichment analysis (GSEA) results for
                            each
                            dataset. In the GSEA tab, the pre-calculated GSEA results are available for users to
                            characterize the functional differences between different cell types. We collected 16,626
                            gene sets from MSigDB, covering KEGG, hallmark. Heatmaps will be shown to display the
                            enriched up-or-down-regulated pathways identified based on differential genes in each
                            cluster.(Beta)</p><br>
                        <img class="selinalogo" src="../assets/doc_gsea.png"><br>

                        <p><b>CCI:</b> In this module, we provide cell-cell interaction(CCI) power by cellchat.(Beta)
                        </p><br>
                        <p><b>Download:</b> You can download the expression matrix and mate information in this
                            module.(Beta)</p><br>
                        <h2 id="annotation">Annotation</h2><br>
                        <p>HUSCH provides cell-type annotation function powered by <a class="doclink"
                                href="https://github.com/SELINA-team/SELINA.py" target="_blank">SELINA</a>, a deep
                            learning-based
                            framework for single cell assignment with multiple reference. You can query data of 54
                            tissue (including Blood, Brain, Lung, Liver, etc) with the three kinds of formats. You can
                            also select <code>Single-cell</code> or <code>Cluster</code> level for user-specific
                            annotation.<br>

                        <p>The annotation result will sent to your email after a while.</p>
                        </p><br>
                        <img class="selinalogo" src="../assets/doc_anno.png"><br>

                        <p>Before using this function, there are a few things to note:</p><br>
                        <p>SELINA supports 3 kinds input format: <code>plain</code>,<code>h5</code> and
                            <code>mtx</code>.
                            You can find
                            detail information in <a class="doclink"
                                href="https://selina.readthedocs.io/en/latest/prepare.html#preprocess-of-query-data"
                                target="_blank">Preprocess-of-query-data</a>
                            part.
                        </p>
                        <p>To guarantee your experience, we only allow uploading <code>.tar.gz</code> files for time
                            saving. You
                            need to compress the original file. Noted that the compressed file should have same prefix
                            with the
                            original file name(eg : <code>tar -zcvf Query.tar.gz Query</code>).</p><br>
                    </div>
                </div>
            </el-tab-pane>
            <el-tab-pane label="Cell type level">
                <el-row :gutter="20">
                    <el-col :span="10">
                        <el-card class="box-card">
                            <el-select v-model="value" placeholder="Select Cell-type Level">
                                <el-option v-for="item in options" :key="item.value" :label="item.label"
                                    :value="item.value" @click.native="selectThing(item)">
                                </el-option>
                            </el-select>
                            <el-scrollbar>
                                <div>
                                    <el-tree ref="tree" accordion :data="treeData" highlight-current
                                        @node-click="getCurrentNode" @node-collapse="collapseNode"
                                        :render-content="renderContent" v-show="isShowlevel2">
                                        <span slot-scope="{ node, data }" class="custom-tree-node">
                                            <span>
                                                <i v-if="data.icon" :class="'tree-icon iconfont ' + data.icon" />
                                                <span class="tree-text">{{ node.label }}</span>
                                            </span>
                                        </span>
                                    </el-tree>
                                </div>
                                <div>
                                    <el-tree ref="tree" accordion :data="treeData3" highlight-current
                                        @node-click="getCurrentNode" @node-collapse="collapseNode"
                                        :render-content="renderContent" v-show="isShowlevel3">
                                        <span slot-scope="{ node, data }" class="custom-tree-node">
                                            <span>
                                                <i v-if="data.icon" :class="'tree-icon iconfont ' + data.icon" />
                                                <span class="tree-text">{{ node.label }}</span>
                                            </span>
                                        </span>
                                    </el-tree>
                                </div>
                            </el-scrollbar>
                        </el-card>
                    </el-col>

                </el-row>
            </el-tab-pane>
        </el-tabs>
    </div>

</template>

<script>
import axios from "axios";
import summary_celltype_l2 from "@/assets/summary_celltype_l2.json"; // 引用
import summary_celltype_l3 from "@/assets/summary_celltype_l3.json"; // 引用

export default {
    data() {
        return {
            treeData: [],
            treeData3: [],
            tempList: [],
            tempList3: [],

            showtissue: false,
            tissuename: "",
            isShowlevel2: false,
            isShowlevel3: false,

            options: [
                {
                    value: "Level2",
                    label: "Major",
                },
                {
                    value: "Level3",
                    label: "Minor",
                },
            ],
            value: "",
        };
    },
    created() {
        this.getTreeData();
        this.getTreeData3();
    },
    methods: {
        selectThing(item) {
            if (item.value === "Level2") {
                //options中的value
                this.isShowlevel2 = true;
                this.isShowlevel3 = false;
            } else if (item.value === "Level3") {
                this.isShowlevel2 = false;
                this.isShowlevel3 = true;
            }
        },

        renderContent(h, { node, data, store }) {
            return (
                <span class="custom-tree-node">
                    <span>{node.label}</span>
                </span>
            );
        },

        getCurrentNode() {
            console.log(this.$refs.tree.getCurrentNode().label);
            this.tissuename = this.$refs.tree.getCurrentNode().label;
            //   console.log(tissuename);
            this.showtissue = true;
        },

        collapseNode() {
            this.showtissue = false;
        },

        async getTreeData() {
            this.tempList = summary_celltype_l2;
            const temp_list = [];

            for (var i = 0; i < this.tempList.length; i++) {
                const obj = {
                    label: this.tempList[i].Tissue,
                    id: i + 1,
                    children: [],
                };
                if (this.tempList[i].celltype.length >= 1) {
                    let temp = this.tempList[i].celltype.length;
                    while (temp > 0) {
                        const obj_child = {
                            label:
                                this.tempList[i].celltype[
                                this.tempList[i].celltype.length - temp
                                ],
                            id: (i + 1) * 100 + this.tempList[i].celltype.length - temp,
                        };
                        temp--;
                        obj.children.push(obj_child);
                    }
                }
                temp_list.push(obj);
            }
            this.treeData = temp_list;
        },

        async getTreeData3() {
            this.tempList3 = summary_celltype_l3;
            const temp_list3 = [];

            for (var i = 0; i < this.tempList3.length; i++) {
                const obj = {
                    label: this.tempList3[i].Tissue,
                    id: i + 1,
                    children: [],
                };
                if (this.tempList3[i].celltype.length >= 1) {
                    let temp = this.tempList3[i].celltype.length;
                    while (temp > 0) {
                        const obj_child = {
                            label:
                                this.tempList3[i].celltype[
                                this.tempList3[i].celltype.length - temp
                                ],
                            id: (i + 1) * 100 + this.tempList3[i].celltype.length - temp,
                        };
                        temp--;
                        obj.children.push(obj_child);
                    }
                }
                temp_list3.push(obj);
            }
            this.treeData3 = temp_list3;
        },
    },
};
</script>
    



<style lang="scss" scoped>
.el-tree {
    background: transparent;
    padding-top: 20px;
}

.custom-tree-node {
    flex: 1;
    display: flex;
    align-items: center;
    justify-content: space-between;
    font-size: 14px;
    padding-right: 8px;
    background-color: blue;
}

.el-tree-node__label {
    font-size: 14px;
}

img {
    width: 50%;
    height: auto;
}

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
</style>