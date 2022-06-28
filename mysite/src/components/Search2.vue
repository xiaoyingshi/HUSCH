<template>
  <div id="search">
    <div id='title'>
      <el-row>
        <el-col>

        </el-col>
      </el-row>
      <h1>Welcome to HUSCH</h1>
      <h3>HUSCH (Human Universal Single Cell Hub) is a scRNA-seq database focusing human universal cell. HUSCH provides detailed cell-type annotation at the single-cell level.Furthermore,we also provide online cell annotation tools and online analysis tools. </h3>
    </div>
    <div id='tab'>
      <el-tabs type="border-card">
      <el-tab-pane label="TIssue Map">
        <el-row>
      <el-card>
                   <div class='test'>
                    <img
                      
                      width="1760px"
                      height="1200px"
                      id='body'
                      class="image"
                      usemap="#ogmap"
                      
                      v-bind:src="
                        require('@/assets/human6.png')
                      "
                    />
                  
                      
                    <map name="ogmap">
                      
                      <area shape="circle" coords="760,112,48" href="#/info_tissue/Brain" alt="Brain">
                      <area shape="circle" coords="880,122,58" href="#/info_tissue/Eye" alt="Eye">
                      <area shape="circle" coords="1004,178,55" href="#/info_tissue/Nose" alt="Nose">
                      <area shape="circle" coords="1116,270,55" href="#/info_tissue/Tonsil" alt="Tonsil">
                      <area shape="circle" coords="1182,390,55" href="#/info_tissue/Lung" alt="Lung">
                      <area shape="circle" coords="1220,510,55" href="#/info_tissue/Heart" alt="Heart">
                      <area shape="circle" coords="1220,630,55" href="#/info_tissue/Spleen" alt="Spleen">
                      <area shape="circle" coords="1212,754,55" href="#/info_tissue/Kideny" alt="Kideny">
                      <area shape="circle" coords="1160,874,55" href="#/info_tissue/Stomach" alt="Stomach">
                      <area shape="circle" coords="1070,986,55" href="#/info_tissue/Intestines" alt="Intestines">
                      <area shape="circle" coords="950,1034,55" href="#/info_tissue/Artery" alt="Artery">
                      <area shape="circle" coords="830,1080,55" href="#/info_tissue/Uterus" alt="Uterus">
                      <area shape="circle" coords="702,1074,55" href="#/info_tissue/" alt="Urinary">
                      <area shape="circle" coords="574,1048,50" href="#/info_tissue/Testis" alt="Testis">
                      <area shape="circle" coords="462,962,50" href="#/info_tissue/Bone-Marrow" alt="Bone-Marrow">
                      <area shape="circle" coords="376,872,50" href="#/info_tissue/Pancreas" alt="Pancreas">
                      <area shape="circle" coords="288,744,55" href="#/info_tissue/Gall-Bladder" alt="Gall-Bladder">
                      <area shape="circle" coords="278,624,55" href="#/info_tissue/Gall-Bladder" alt="Liver">
                      <area shape="circle" coords="286,500,55" href="#/info_tissue/Breast" alt="Breast">
                      <area shape="circle" coords="326,378,55" href="#/info_tissue/Muscle" alt="Muscle">
                      <area shape="circle" coords="404,266,55" href="#/info_tissue/" alt="Thymus- Gland">
                      <area shape="circle" coords="512,168,55" href="#/info_tissue/Thyroid" alt="Thyroid">
                      <area shape="circle" coords="640,128,55" href="#/info_tissue/Esophagus" alt="Esophagus">
                      <area shape="circle" coords="1530,64,50" href="#/info_tissue/Blood" alt="Blood">
                      <area shape="circle" coords="1530,206,50" href="#/info_tissue/Skin" alt="Skin">
                      <area shape="circle" coords="1530,360,50" href="#/info_tissue/Airway" alt="Airway">
                      <area shape="circle" coords="1530,520,50" href="#/info_tissue/Adipose" alt="Adipose">
                      <area shape="circle" coords="1530,670,50" href="#/info_tissue/Adrenal-Gland" alt="Adrenal-Gland">
                      <area shape="circle" coords="1530,820,50" href="#/info_tissue/Bladder" alt="Bladder">
                      <area shape="circle" coords="1530,968,50" href="#/info_tissue/Placenta" alt="Placenta">
                      <area shape="circle" coords="1530,1108,50" href="#/info_tissue/Falloplan" alt="Falloplan">


                    </map>
                    <!-- </el-col>
        </el-row>
      </el-card> -->
      </div>
      </el-card>
    </el-row>
      </el-tab-pane>
    <el-tab-pane label="Tissue List">
      <div id='list'>
      <el-row>
        <el-card>
        <el-divider content-position="left"><span style="color: rgb(26, 96, 67) ;">Tissue List</span></el-divider>
        <el-col :span = '24'>
                
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 2</p>
                      <p>CellNumber:</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Eosinophils,Mast,Mono/Macro,Myofibroblasts,Natural killer cell(NK),Neutrophils,Plasma,Regulatory T cell(Treg),TMKI67,pDC</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Adrenal-Gland"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Adrenal-Gland.png" alt="">Adrenal-Gland</div></router-link>
                    </el-popover>
                    <!-- <router-link to="info_tissue/Adrenal-Gland"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Adrenal-Gland.png" alt="">Adrenal-Gland</div></router-link> -->
                  </el-col>
                  
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 2</p>
                      <p>CellNumber:</p>
                      <p>Celltype : </p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Airway"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Airway.png" alt="">Airway</div></router-link>
                    </el-popover>
                    
                  </el-col>
                  <el-col :span = '5'>
                    
                    
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 2</p>
                      <p>CellNumber:</p>
                      <p>Celltype :  B,Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Myofibroblast,T</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Artery"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Artery.png" alt="">Artery</div></router-link>
                    </el-popover>
                    
                  </el-col>
                  
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 5</p>
                      <p>CellNumber:</p>
                      <p>Celltype :  B,Basal,CD4T,CD8T,Endothelial,Epithelial,Fibroblast,Intermediate(Umbrella-Basal),Interstitial,Macrophage,Mast,Monocyte,Muscle,Natural killer cell(NK),Neutrophil,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Bladder"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Bladder.png" alt="">Bladder</div></router-link>
                    </el-popover>
                    
                  </el-col>
                

              </el-col>
              <br>
        <el-col :span = '24'>
                
                  
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 13</p>
                      <p>CellNumber:</p>
                      <p>Celltype :  B,B/NK Progenitor,Basophil,CD4T,CD8T,Common myeloid progenitor(CMP),Dendritic cell(DC),EPC,Erythroblast,Erythroid,Granulocyte_monocyte progenitor(GMP),Hematopoietic stem cell(HSC),Hemopoietic stem/progenitor cell(HSPC),Lymphoid Progenitor,Lymphoid-primed multipotent progenitor(LMPP),MPP,Macrophage,Megakaryocyte-erythroid progenitor(MEP),Monocyte,Multi-lymphoid progenitor(MLP),Natural killer cell(NK),Neutrophil,T Progenitor</p>
                      
                    </div>
                    
                    <router-link slot="reference" to="info_tissue/Bone-Marrow"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Bone-Marrow.png" alt="">Bone-Marrow</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 15</p>
                      <p>CellNumber:</p>
                      <p>Celltype : </p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Blood"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Blood.png" alt="">Blood</div></router-link>
                    </el-popover>
                    
                  </el-col>
                  
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 9</p>
                      <p>CellNumber:</p>
                      <p>Celltype :  B,CD4T,CD8T,Endothelial,Epithelial,Excitatory Neuron,Fibroblast,Glial,Inhibitory Neuron,Interneuron,Macrophage,Mast,Muscle,Natural killer cell(NK),Neutrophil,Pericyte,Plasma</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Brain"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Brain.jpg" alt="">Brain</div></router-link>
                    </el-popover>
                    
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 5</p>
                      <p>CellNumber:</p>
                      <p>Celltype :  Basal,Luminal</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Breast"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Breast.png" alt="">Breast</div></router-link>
                    </el-popover>
                    
                  </el-col>
                

        </el-col>
        <br>
        <el-col :span = '24'>
                
                  
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber:</p>
                      <p>Celltype :   B,CD4T,CD8T,Endothelial,Fibroblast,Macrophage,Melanocyte,Pericyte,Photoreceptor/RPE,Red Blood Cell(RBC),Schwann,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Choroid"><div id="tiss" type="primary" round>Choroid</div></router-link>
                    </el-popover>
                    
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber:</p>
                      <p>Celltype :  B,CD4T,CD8T,Endothelial,Fibroblast,Macrophage,Melanocyte,Pericyte,Photoreceptor/RPE,Red Blood Cell(RBC),Schwann,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Esophagus"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Esophagus.png" alt="">Esophagus</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 3</p>
                      <p>CellNumber:</p>
                      <p>Celltype : Amacrine,B,Bipolar,CD4T,CD8T,Cone,Endothelial,Epithelial,Fibroblast,Ganglion,Glial,Horizontal,Macrophage,Mast,Melanocyte,Natural killer cell(NK),Pericyte,Rod</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Eye"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Eye.png" alt="">Eye</div></router-link>
                     </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber:</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Neutrophil</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Fallopian"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Fallopian.png" alt="">Fallopian</div></router-link>
                    </el-popover>
                  </el-col>
                

        </el-col>
        <br>
        <el-col :span = '24'>
                
                  
                  <el-col :span = '5'>
                    
                    
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 2</p>
                      <p>CellNumber:</p>
                      <p>Celltype :B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Natural killer cell(NK),Neutrophil</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Gall-Bladder"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Gall-Bladder.png" alt="">Gall-Bladder</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                     <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 2</p>
                      <p>CellNumber:</p>
                      <p>Celltype :Dendritic cell(DC),Endothelial,Fibroblast,Macrophage,Mast,Muscle,Neutrophil,T</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Heart"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Heart.jpg" alt="">Heart</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 14</p>
                      <p>CellNumber:</p>
                      <p>Celltype :B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Glial,Innate Lymphoid Cell(ILC),Intraepithelial Lymphocyte(IEL),Macrophage,Mast,Mesenchymal,Monocyte,Muscle,Myofibroblast,Natural killer cell(NK),Neuron,Neutrophil,Paneth-like,Pericyte,Plasma,Stem,Transit amplifying(TA)</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Intestine"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Intestines.jpg" alt="">Intestine</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 14</p>
                      <p>CellNumber:</p>
                      <p>Celltype :B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Monocyte,Mononuclear Phagocyte(MNP),Muscle,Myofibroblast,NKT,Natural killer cell(NK),Neutrophil,Pericyte,Plasma,Regulatory T cell(Treg),Stem,Tfh,Tprolif</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Kidney"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Kidney.jpg" alt="">Kidney</div></router-link>
                    </el-popover>
                  </el-col>
                

        </el-col>
        <br>
        <el-col :span = '24'>
                
                  
                  <el-col :span = '5'>
                    
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 2</p>
                      <p>CellNumber:</p>
                      <p>Celltype :B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Erthyroid,Kupffer,Macrophage,Mast,Mesenchymal,Monocyte,Muscle,Natural killer cell(NK),Neutrophil,Plasma</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Liver"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Liver.jpg" alt="">Liver</div></router-link>
                    </el-popover>
                    
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 2</p>
                      <p>CellNumber:</p>
                      <p>Celltype :B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Erthyroid,Kupffer,Macrophage,Mast,Mesenchymal,Monocyte,Muscle,Natural killer cell(NK),Neutrophil,Plasma</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Lung"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Lung.jpg" alt="">Lung</div></router-link>
                  </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                     <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 2</p>
                      <p>CellNumber:</p>
                      <p>Celltype :Adipocyte,B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Fibroblast,Macrophage,Mast,Muscle Progenitor,Natural killer cell(NK),Neutrophil,Pericyte,Platelet,Skeletal Muscle,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Muscle"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Muscle.png" alt="">Muscle</div></router-link>
                   </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber:</p>
                      <p>Celltype :B,CD4T,CD8T,Dendritic cell(DC),Epithelial,Fibroblast,Macrophage,Mast,Monocyte,Natural killer cell(NK),Neuron,Olfactory ensheathing Glia(OEC),Pericyte,Plasma,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Nose"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Nose.png" alt="">Nose</div></router-link>
                    </el-popover>
                  </el-col>
                

        </el-col>
        <br>
        <el-col :span = '24'>
                
                  
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 4</p>
                      <p>CellNumber:</p>
                      <p>Celltype :Endocrine,Endothelial,Epithelial,Exocrine,Fibroblast,Macrophage,Mast,Mesenchymal,Muscle,Pancreatic Stellate,Schwann</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Pancreas"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Pancreas.jpg" alt="">Pancreas</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber:</p>
                      <p>Celltype :B,CD4T,CD8T,Endothelial,Fibroblast,Mesothelial,Mononuclear Phagocyte,Myofibroblast</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Lung"><div id="tiss" type="primary" round>Peritoneum</div></router-link>
                  </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 5</p>
                      <p>CellNumber:</p>
                      <p>Celltype :Antigen-presenting,CD4T,CD8T,Decidual,Dendritic cell(DC),Endothelial,Epithelial,Erythroid,Fibroblast,Granulocyte,Hofbauer,Macrophage,Monocyte,Muscle,Neutrophil,Proliferating,Trophoblast,Vascular Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Placenta"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Placenta.png" alt="">Placenta</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : </p>
                      <p>CellNumber:</p>
                      <p>Celltype :Epithelial,Fibroblast,Macrophage,Natural killer cell(NK),Neutrophil</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Plerua"><div id="tiss" type="primary" round>Plerua</div></router-link>
                   </el-popover>
                  </el-col>
                

        </el-col>
        <br>
        <el-col :span = '24'>
                
                  
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 3</p>
                      <p>CellNumber:</p>
                      <p>Celltype :CD4T,CD8T,Endothelial,Epithelial,Fibroblast,Leukocyte,Macrophage,Muscle,Natural killer cell(NK),Neuroendocrine,Neutrophil,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Prostate"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Urinary.jpg" alt="">Prostate</div></router-link>
                     </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 2</p>
                      <p>CellNumber:</p>
                      <p>Celltype :Bulge,CD4T,CD8T,Endothelial,Epidermal Stem,Fibroblast,IFE Basal,IFE Granular,IFE Mitotic,IFE Spinous,IRS/cuticle,Infundibulum,Isthmus,Keratinocyte,Langerhans,Macrophages/DC,Matrix/cortex/medulla,Melanocyte,Mesenchymal,ORS B,Pericyte,Sebaceous/apocrine</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Skin"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Skin.png" alt="">Skin</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber:</p>
                      <p>Celltype :B,CD4T,CD8T,Endothelial,Erythroid,Lymphoid Progenitor,Macrophage,Neutrophil</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Spleen"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Spleen.png" alt="">Spleen</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 3</p>
                      <p>CellNumber:</p>
                      <p>Celltype :B,CD4T,CD8T,Endothelial,Epithelial,Fibroblast,Granulocyte,Macrophage,Mast,Muscle,Natural killer cell(NK),Neuroendocrine</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Stomach"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Stomach.jpg" alt="">Stomach</div></router-link>
                    </el-popover>
                  </el-col>
                  
                

        </el-col>
        <br>  
        <el-col :span = '24'>
                
                  
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 4</p>
                      <p>CellNumber:</p>
                      <p>Celltype :Apical/Narrow,Basal,Clear,Efferent Duct,Endothelial,Germ,Leydig,Macrophage,Peritubular Myoid Cell(PTM),Principal,Sertoli,Sperm,Spermatid,Spermatocyte,Spermatogonia,Spermatogonial Stem cell(SSC),Stromal,Vascular Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Testis"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Testis.jpg" alt="">Testis</div></router-link>
                     </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 2</p>
                      <p>CellNumber:</p>
                      <p>Celltype :B,CD4T,CD8T,Endothelial,Epithelial,Fibroblast,Granulocyte,Macrophage,Mast,Muscle,Natural killer cell(NK),Neuroendocrine</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Thyroid"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Thyroid.jpg" alt="">Thyroid</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber:</p>
                      <p>Celltype :Dendritic cell(DC)e</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Tonsil"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Tonsil.png" alt="">Tonsil</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber:</p>
                      <p>Celltype :B,CD4T,Endometrial,Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Natural killer cell(NK),Regulatory T cell(Treg),Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Ureter"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Urinary.jpg" alt="">Ureter</div></router-link>
                    </el-popover>
                  </el-col>
                  
                

        </el-col>
        <br>
        <el-col :span = '24'>
                
                  
                  
                  <el-col :span = '5'>
                    <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 5</p>
                      <p>CellNumber:</p>
                      <p>Celltype :Adipocyte,B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Granulocyte,Macrophage,Mast,Monocyte,Natural killer cell(NK),Neutrophil,Stem</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Visceral-Adipose"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Adipose.png" alt="">Visceral-Adipose</div></router-link>
                    </el-popover>
                  </el-col>
                 

        </el-col>
      </el-card>
      </el-row>
      </div>
     
    </el-tab-pane>
    </el-tabs>
    </div>
    
    
    
       
      
     
              
    <!-- </el-row> -->

  </div>
</template>


<script>
import celltype_json from "@/assets/summary_celltype_l2.json"; // 引用
export default {
  data() {
    return {
      // img_src :'@/assets/human6.png?t='+Math.random()
      celltype_list:celltype_json
    }
  },
  created() {
    console.log('creat')
    // location.reload();
  },
  activated() {
    // location.reload();
    console.log('acti')
    
    // if (location.href.indexOf("#reloaded") == -1) {
    //     location.href = location.href + "#reloaded";
    //     location.reload();
    // }
  },
  methods: {
    submitUpload() {
      this.$refs.upload.submit();
    },
  },
};
</script>
<style>
.el-row{
  width: 1760px;
}
.el-col{
  margin: 5px;
}
.el-buttom{
  background-color: rgb(26, 96, 67);
}
.selinalogo {
    height: 200px;
    padding: 10px;
    /* padding-right: 0px;
    padding-bottom: 10px;
    padding-top: 3px; */
  }
#tiss{
  display: flex;
  flex-direction: row;
  justify-content:flex-start;
  align-items:center;
  border-radius: 30px;
  /* border: 3px solid rgb(26, 96, 67); */
  padding: 10px;
  height: 70px;
  font-size: 40px;
  color: rgb(33, 56, 50);
  background-color:white;
  
}

/* #body{
  width:"1760px" !important;
  height:"1200px";
} */
.el-card{
  width:"1760px";
}
#search{
  display: flex;
  flex-direction: column;
  justify-content: center;
  /* width: 2010px; */
  align-items:center;

}
.el-divider{
    color: rgb(26, 96, 67);
    background-color: rgb(26, 96, 67);
    font-size: 70px;
}
/* img{
  width:"1760px";
  height:"1200px";
  object-fit: none;
} */
#title{
  display: flex;
  flex-direction: column;
  justify-content: center;
  align-items:center;
  color: rgb(26, 96, 67);
  width: 1500px;
  margin-bottom: 80px;
}
.el-tabs{
  width: 1790px;
}

</style>


