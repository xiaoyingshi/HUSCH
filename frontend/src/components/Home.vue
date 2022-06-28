<template>
  <div id="home">
    <el-container class="home-container">
      <el-header height="70px">
        <!-- <div class="ico"> -->

        <img
          class="selinalogo"
          src="../assets/husch5.png"
          style="height: 70px; padding: 5px"
        />
        <p style="font-size: 40px">Human Universal Single Cell Hub</p>
        <!-- </div> -->
        <!-- <div>
          <i class="el-icon-message"></i
          ><span style="font-size: 15px; color: white"
            >bioinfo.sxy@gmail.com</span
          >
        </div> -->
      </el-header>

      <el-container>
        <el-aside :width="isCollapse ? '64px' : '160px'">
          <el-scrollbar style="height: 100%">
            <div class="toggle-button" @click="togleCollapse">|||</div>
            <el-menu
              unique-opened
              :collapse="isCollapse"
              :collapse-transition="false"
              default-active="1"
              background-color="#dce1de"
              text-color="#567568"
              router
            >
              <el-menu-item
                index="/home"
                @mouseenter="changeActive($event)"
                @mouseleave="removeActive($event)"
              >
                <i class="el-icon-house"></i>
                <span slot="title">Home</span>
              </el-menu-item>
              <el-menu-item index="/info">
                <i class="el-icon-shopping-cart-2"></i>
                <span slot="title">Datasets</span>
              </el-menu-item>
              <el-menu-item index="/gene">
                <i class="el-icon-orange"></i>
                <span slot="title">Gene</span>
              </el-menu-item>
              <el-menu-item index="/annotation">
                <i class="el-icon-search"></i>
                <span slot="title">Annotation</span>
              </el-menu-item>
              <!-- <el-menu-item index="/Statistics">
                <i class="el-icon-data-line"></i>
                <span slot="title">Statistics</span>
              </el-menu-item> -->
              <el-menu-item index="/documentation">
                <i class="el-icon-document"></i>
                <span slot="title">Documentation</span>
              </el-menu-item>
            </el-menu>
          </el-scrollbar>
        </el-aside>
        <el-container>
          <el-main>
            <router-view></router-view>
            <!-- <keep-alive> <router-view></router-view></keep-alive> -->
          </el-main>
          <el-footer></el-footer>
        </el-container>
      </el-container>
    </el-container>
  </div>
</template>



<script>
export default {
  data() {
    return {
      isCollapse: false,
    };
  },
  methods: {
    togleCollapse() {
      this.isCollapse = !this.isCollapse;
    },
    changeActive($event) {
      $event.currentTarget.className = "animate__animated animate__headShake";
    },
    removeActive($event) {
      $event.currentTarget.className = "";
    },
  },
};
</script>

<style  lang="scss" scoped>
.selinalogo {
  height: 50px;
  padding: 10px;
  /* padding-right: 0px;
    padding-bottom: 10px;
    padding-top: 3px; */
}
.home-container {
  height: 100vh;
}
.ico {
  display: flex;
  flex-direction: row;
}
.el-menu-item {
  font-weight: bold;
  color: #567568;

}


.el-header {
  background-color: #dce1de;
  display: flex;
  flex-wrap: wrap;
  padding-left: 0;
  align-items: center;
  color: #567568;
  font-weight: bold;
  font-family: Arial;
  border-bottom: 2px solid #dce1de;
}


.el-header > p {
  margin-left: 80px;
}

.el-footer {
  background-color: #dce1de;
}

.el-aside {
  background-color: #dce1de;
}
.el-aside.el-menu {
  border: none;
  font-weight: bold;
  height: 100%;
}

.el-main {
  background-color: white;
}

.logo {
  width: 70px;
}

.toggle-button {
  background-color: #567568;
  font-size: 10px;
  line-height: 22px;
  font-weight: bold;
  color: #dce1de;
  text-align: center;
  letter-spacing: 0.2em;

  cursor: pointer;
}

.el-icon-message {
  font-size: 20px;
  font-weight: bold;
  color: #567568;
  padding-right: -2px;
}
</style>