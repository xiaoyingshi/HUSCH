<template>
  <div id="search">
    <div id='title'>
      
      
      <h1>Welcome to HUSCH</h1>
      <h3>HUSCH (Human Universal Single Cell Hub) is a scRNA-seq database focusing human universal cell. HUSCH provides detailed cell-type annotation at the single-cell level.Furthermore,we also provide online cell annotation tools and online analysis tools. </h3>
    </div>

    

 
 <div id='cont'>
  <el-card>
        <el-divider content-position="left"><span style="color: rgb(26, 96, 67) ;">Tissue List</span></el-divider>
       <!-- <div id='tissue'> -->
         <el-row>

        
      
      
          <el-col :span = '6'>
            <div class="animate__animated animate__backInLeft">
                <img src="@/assets/human7.png" alt="">
            </div>
           
          </el-col>
          
       
         <el-col :span = '14'>
                 
        <el-row >
                <el-col :span = '5'>
                    
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 4</p>
                      <p>CellNumber : 16593</p>
                      <p>Celltype : Adipocyte,B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Granulocyte,Macrophage,Mast,Monocyte,Natural killer cell(NK),Neutrophil,Stem</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Adipose"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Adipose.png" alt="">Adipose</div></router-link>
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
                      <p>CellNumber: 16197</p>
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
                      <p>Dataset : 3</p>
                      <p>CellNumber : 13647</p>
                      <p>Celltype : Ciliated,Basal,Secretory,Brush+PNEC</p>
                      
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
                      <p>Dataset : 1</p>
                      <p>CellNumber : 9652</p>
                      <p>Celltype :  B,Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Myofibroblast,T</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Artery"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Artery.png" alt="">Artery</div></router-link>
                    </el-popover>
                    
                  </el-col>
                  
                  
                

        <!-- </el-row> -->
              <br>
        <!-- <el-row> -->
                <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 4</p>
                      <p>CellNumber : 17087</p>
                      <p>Celltype :  B,Basal,CD4T,CD8T,Endothelial,Epithelial,Fibroblast,Intermediate(Umbrella-Basal),Interstitial,Macrophage,Mast,Monocyte,Muscle,Natural killer cell(NK),Neutrophil,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Bladder"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Bladder.png" alt="">Bladder</div></router-link>
                    </el-popover>
                    
                  </el-col>
                   <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 13</p>
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
                      <p>Dataset : 1</p>
                      <p>CellNumber:</p>
                      <p>Celltype :  B,B/NK Progenitor,Basophil,CD4T,CD8T,Common myeloid progenitor(CMP),Dendritic cell(DC),EPC,Erythroblast,Erythroid,Granulocyte_monocyte progenitor(GMP),Hematopoietic stem cell(HSC),Hemopoietic stem/progenitor cell(HSPC),Lymphoid Progenitor,Lymphoid-primed multipotent progenitor(LMPP),MPP,Macrophage,Megakaryocyte-erythroid progenitor(MEP),Monocyte,Multi-lymphoid progenitor(MLP),Natural killer cell(NK),Neutrophil,T Progenitor</p>
                      
                    </div>
                    
                    <router-link slot="reference" to="info_tissue/Bone"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Bone.png" alt="">Bone</div></router-link>
                    </el-popover>
                  </el-col>
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
                    
                    <router-link slot="reference" to="info_tissue/Bone-Marrow"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Bone-Marrow.png" alt="">Bone Marrow</div></router-link>
                    </el-popover>
                  </el-col>
                 
                  
                  
                

        <!-- </el-row> -->
        <br>
        <!-- <el-row> -->
                <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 9</p>
                      <p>CellNumber : 88228</p>
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
                      <p>Dataset : 4</p>
                      <p>CellNumber : 25685</p>
                      <p>Celltype :  Basal,Luminal</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Breast"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Breast.png" alt="">Breast</div></router-link>
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
                      <p>CellNumber : 4587</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Neutrophil</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Cervix"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Cervix.png" alt="">Cervix</div></router-link>
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
                      <p>CellNumber : 9839</p>
                      <p>Celltype :   B,CD4T,CD8T,Endothelial,Fibroblast,Macrophage,Melanocyte,Pericyte,Photoreceptor/RPE,Red Blood Cell(RBC),Schwann,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Choroid"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Eye.png" alt="">Choroid</div></router-link>
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
                      <p>CellNumber : 11024</p>
                      <p>Celltype :  B,CD4T,CD8T,Endothelial,Fibroblast,Macrophage,Melanocyte,Pericyte,Photoreceptor/RPE,Red Blood Cell(RBC),Schwann,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Common.bile.duct"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Common.bile.duct.png" alt="">Common Bile Duct</div></router-link>
                    </el-popover>
                  </el-col>
                  
                

        <!-- </el-row> -->
        <br>
        <!-- <el-row> -->
          <el-col :span = '5'>
                   <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber : 4587</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Neutrophil</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Decidua"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Decidua.png" alt="">Decidua</div></router-link>
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
                      <p>CellNumber : 11024</p>
                      <p>Celltype :  B,CD4T,CD8T,Endothelial,Fibroblast,Macrophage,Melanocyte,Pericyte,Photoreceptor/RPE,Red Blood Cell(RBC),Schwann,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Ear"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Ear.png" alt="">Ear</div></router-link>
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
                      <p>CellNumber : 11024</p>
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
                      <p>Dataset : 6</p>
                      <p>CellNumber : 87008</p>
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
                      <p>CellNumber : 4587</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Neutrophil</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Fallopian"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Fallopian.png" alt="">Fallopian</div></router-link>
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
                      <p>CellNumber : 13936</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Natural killer cell(NK),Neutrophil</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Gall-Bladder"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Gall-Bladder.png" alt="">Gall-Bladder</div></router-link>
                    </el-popover>
                  </el-col>
                  
                

        <!-- </el-row> -->
        <br>
        <!-- <el-row> -->
                
                  <el-col :span = '5'>
                     <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 3</p>
                      <p>CellNumber : 2786</p>
                      <p>Celltype : Dendritic cell(DC),Endothelial,Fibroblast,Macrophage,Mast,Muscle,Neutrophil,T</p>
                      
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
                      <p>Dataset : 11</p>
                      <p>CellNumber : 112727</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Monocyte,Mononuclear Phagocyte(MNP),Muscle,Myofibroblast,NKT,Natural killer cell(NK),Neutrophil,Pericyte,Plasma,Regulatory T cell(Treg),Stem,Tfh,Tprolif</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Kidney"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Kidney.jpg" alt="">Kidney</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 7</p>
                      <p>CellNumber : 44812</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Erthyroid,Kupffer,Macrophage,Mast,Mesenchymal,Monocyte,Muscle,Natural killer cell(NK),Neutrophil,Plasma</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Liver"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Liver.jpg" alt="">Liver</div></router-link>
                    </el-popover>
                    
                  </el-col>
                  
                

        <!-- </el-row> -->
        <br>
        <!-- <el-row> -->
                
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 12</p>
                      <p>CellNumber : 114267</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Erthyroid,Kupffer,Macrophage,Mast,Mesenchymal,Monocyte,Muscle,Natural killer cell(NK),Neutrophil,Plasma</p>
                      
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
                      <p>Dataset : 1</p>
                      <p>CellNumber : 29428</p>
                      <p>Celltype : Adipocyte,B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Fibroblast,Macrophage,Mast,Muscle Progenitor,Natural killer cell(NK),Neutrophil,Pericyte,Platelet,Skeletal Muscle,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Lymph-Node"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Lymph-Node.png" alt="">Lymph Node</div></router-link>
                   </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                     <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 4</p>
                      <p>CellNumber : 29428</p>
                      <p>Celltype : Adipocyte,B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Fibroblast,Macrophage,Mast,Muscle Progenitor,Natural killer cell(NK),Neutrophil,Pericyte,Platelet,Skeletal Muscle,Smooth Muscle</p>
                      
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
                      <p>Dataset : 2</p>
                      <p>CellNumber : 33060</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Epithelial,Fibroblast,Macrophage,Mast,Monocyte,Natural killer cell(NK),Neuron,Olfactory ensheathing Glia(OEC),Pericyte,Plasma,Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Nose"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Nose.png" alt="">Nose</div></router-link>
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
                      <p>CellNumber : 4587</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Neutrophil</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Ovary"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Ovary.png" alt="">Ovary</div></router-link>
                    </el-popover>
                  </el-col> 
                

        <!-- </el-row> -->
        <br>
        <!-- <el-row> -->
                
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 4</p>
                      <p>CellNumber : 21845</p>
                      <p>Celltype : Endocrine,Endothelial,Epithelial,Exocrine,Fibroblast,Macrophage,Mast,Mesenchymal,Muscle,Pancreatic Stellate,Schwann</p>
                      
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
                      <p>CellNumber : 17198</p>
                      <p>Celltype : B,CD4T,CD8T,Endothelial,Fibroblast,Mesothelial,Mononuclear Phagocyte,Myofibroblast</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Peritoneum"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Peritoneum.png" alt="">Peritoneum</div></router-link>
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
                      <p>CellNumber : 65500</p>
                      <p>Celltype : Antigen-presenting,CD4T,CD8T,Decidual,Dendritic cell(DC),Endothelial,Epithelial,Erythroid,Fibroblast,Granulocyte,Hofbauer,Macrophage,Monocyte,Muscle,Neutrophil,Proliferating,Trophoblast,Vascular Smooth Muscle</p>
                      
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
                      <p>Dataset : 1</p>
                      <p>CellNumber : 9942</p>
                      <p>Celltype : Epithelial,Fibroblast,Macrophage,Natural killer cell(NK),Neutrophil</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Pleura"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Pleura.png" alt="">Pleura</div></router-link>
                   </el-popover>
                  </el-col>
                  
                  
                

        <!-- </el-row> -->
        <br>  
        <!-- <el-row > -->
                <!-- <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber : 25247</p>
                      <p>Celltype : Bulge,CD4T,CD8T,Endothelial,Epidermal Stem,Fibroblast,IFE Basal,IFE Granular,IFE Mitotic,IFE Spinous,IRS/cuticle,Infundibulum,Isthmus,Keratinocyte,Langerhans,Macrophages/DC,Matrix/cortex/medulla,Melanocyte,Mesenchymal,ORS B,Pericyte,Sebaceous/apocrine</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Rectum"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Rectum.png" alt="">Rectum</div></router-link>
                    </el-popover>
                  </el-col> -->
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 3</p>
                      <p>CellNumber : 15321</p>
                      <p>Celltype : CD4T,CD8T,Endothelial,Epithelial,Fibroblast,Leukocyte,Macrophage,Muscle,Natural killer cell(NK),Neuroendocrine,Neutrophil,Smooth Muscle</p>
                      
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
                      <p>Dataset : 11</p>
                      <p>CellNumber : 25247</p>
                      <p>Celltype : Bulge,CD4T,CD8T,Endothelial,Epidermal Stem,Fibroblast,IFE Basal,IFE Granular,IFE Mitotic,IFE Spinous,IRS/cuticle,Infundibulum,Isthmus,Keratinocyte,Langerhans,Macrophages/DC,Matrix/cortex/medulla,Melanocyte,Mesenchymal,ORS B,Pericyte,Sebaceous/apocrine</p>
                      
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
                      <p>Dataset : 2</p>
                      <p>CellNumber : 15664</p>
                      <p>Celltype : B,CD4T,CD8T,Endothelial,Erythroid,Lymphoid Progenitor,Macrophage,Neutrophil</p>
                      
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
                      <p>Dataset : 4</p>
                      <p>CellNumber : 12568</p>
                      <p>Celltype : B,CD4T,CD8T,Endothelial,Epithelial,Fibroblast,Granulocyte,Macrophage,Mast,Muscle,Natural killer cell(NK),Neuroendocrine</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Stomach"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Stomach.jpg" alt="">Stomach</div></router-link>
                    </el-popover>
                  </el-col>
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 16</p>
                      <p>CellNumber : 119064</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Glial,Innate Lymphoid Cell(ILC),Intraepithelial Lymphocyte(IEL),Macrophage,Mast,Mesenchymal,Monocyte,Muscle,Myofibroblast,Natural killer cell(NK),Neuron,Neutrophil,Paneth-like,Pericyte,Plasma,Stem,Transit amplifying(TA)</p>
                      
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
                      <p>Dataset : 1</p>
                      <p>CellNumber : 59857</p>
                      <p>Celltype : Apical/Narrow,Basal,Clear,Efferent Duct,Endothelial,Germ,Leydig,Macrophage,Peritubular Myoid Cell(PTM),Principal,Sertoli,Sperm,Spermatid,Spermatocyte,Spermatogonia,Spermatogonial Stem cell(SSC),Stromal,Vascular Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Thymus"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Thymus.jpg" alt="">Thymus</div></router-link>
                     </el-popover>
                  </el-col>
                

        <!-- </el-row> -->
        <br>
        <!-- <el-row> -->
                
                
                  <!-- <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber : 59857</p>
                      <p>Celltype : Apical/Narrow,Basal,Clear,Efferent Duct,Endothelial,Germ,Leydig,Macrophage,Peritubular Myoid Cell(PTM),Principal,Sertoli,Sperm,Spermatid,Spermatocyte,Spermatogonia,Spermatogonial Stem cell(SSC),Stromal,Vascular Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Thymus"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Thymus.jpg" alt="">Thymus</div></router-link>
                     </el-popover>
                  </el-col> -->
                  <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 2</p>
                      <p>CellNumber : 12486</p>
                      <p>Celltype : B,CD4T,CD8T,Endothelial,Epithelial,Fibroblast,Granulocyte,Macrophage,Mast,Muscle,Natural killer cell(NK),Neuroendocrine</p>
                      
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
                      <p>CellNumber : 6000</p>
                      <p>Celltype :Dendritic cell(DC)e</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Tonsil"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Tonsil.png" alt="">Tonsil</div></router-link>
                    </el-popover>
                  </el-col>
                  <!-- </el-row> -->
                  
                  <!-- <el-row> -->

                  
                  <el-col :span = '5'>
                    
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber : 1886</p>
                      <p>Celltype :B,CD4T,Endometrial,Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Natural killer cell(NK),Regulatory T cell(Treg),Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Ureter"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Urinary.jpg" alt="">Ureter</div></router-link>
                    </el-popover>
                  </el-col>
                  
                   <!-- </el-row> -->
                  <!-- <el-row> -->

                    <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 4</p>
                      <p>CellNumber : 59857</p>
                      <p>Celltype : Apical/Narrow,Basal,Clear,Efferent Duct,Endothelial,Germ,Leydig,Macrophage,Peritubular Myoid Cell(PTM),Principal,Sertoli,Sperm,Spermatid,Spermatocyte,Spermatogonia,Spermatogonial Stem cell(SSC),Stromal,Vascular Smooth Muscle</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Testis"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Testis.jpg" alt="">Testis</div></router-link>
                     </el-popover>
                  </el-col>
                    <!-- </el-row> -->
        <!-- <el-row> -->

                  
                 
               <el-col :span = '5'>
                    <el-popover
                    placement="right"
                    width="400"
                    trigger="hover"
                    >
                    <div>
                      <p>Dataset : 1</p>
                      <p>CellNumber : 4587</p>
                      <p>Celltype : B,CD4T,CD8T,Dendritic cell(DC),Endothelial,Epithelial,Fibroblast,Macrophage,Mast,Muscle,Neutrophil</p>
                      
                    </div>
                    <router-link slot="reference" to="info_tissue/Uterus"><div id="tiss" type="primary" round><img height="70px" width="70px" src="@/assets/Uterus.png" alt="">Uterus</div></router-link>
                    </el-popover>
                  </el-col> 
                 
                 
                   
                  
                  
                  
                  
                 
                  
        </el-row>
         
          </el-col>
      </el-card>
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
<style scoped>

.el-row{
  width: 120%;
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



#search{
  display: flex;
  flex-direction: column;
  justify-content: center;
  
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
  /* color: rgb(55,123,213); */
  width: 70%;
  margin-bottom: 80px;
}

#tissue{
  /* display: flex;
  flex-direction: row; */
  /* justify-content: space-between; */
  
}
#tiss{
  display: flex;
  flex-direction: row;
  align-items:center;
  font-size: 35px;
  height: 70px;
  color: rgb(33, 56, 50);
}
.el-card{
  width: 100%;
}
#cont{
  width:100%;
}
#body{
  width: 20%;
}
#lisr{
  width: 80%;
}

.link {
    color: black;
    font-weight: bold;
}

/* CSS link color (red) */
.link:hover {
    color: rgb(3, 138, 81);
}


</style>


