<template>
  <div id="app" >
    <keep-alive>
      <router-view :key="this.$route.path"/>
          <router-link to="/user">User</router-link>
          Error
    </keep-alive>
    <!-- <Home/> -->
  </div>
</template>

<script>
import Home from "./components/Home.vue";
import _ from 'lodash'
export default {
  name: "app",
  components: {
    Home,
  },
};
//  var app = new Vue({
//    el: '#app',
//    components: {
//       Home
//       }
//    })
</script>


<style>
@import "./styles/base.css";
</style>
