import '@babel/polyfill'
import 'mutationobserver-shim'
import Vue from 'vue'
import App from './App.vue'
import router from './router'
import axios from 'axios'
import VueAxios from 'vue-axios'
import animate from "animate.css"
Vue.use(animate)


// import VueGtag from 'vue-gtag'

// import VueResource from 'vue-resource'
// import Axios from 'axios'
// import "./assets/icon/iconfont.css"

// require('./mock');
Vue.use(VueAxios, axios)
// Vue.use(VueResource)
Vue.prototype.$axios = Axios

import './plugins/element.js'
// import ELEMENT from 'element-ui';
// Vue.use(ELEMENT, {size: 'small', zIndex: 3000});


import './plugins/bootstrap-vue'

import qs from 'qs';  
import { Axios } from 'axios'
Vue.prototype.$qs = qs;

new Vue({
  router, // 挂载到Vue实例
  render: h => h(App),
}).$mount('#app')



Vue.config.productionTip = false
// Vue.use(vueAxios,axios)
// Vue.use(VueGtag, {
//   config: { id: "husch-344811" } 
// }, router);