import Vue from 'vue'
import Router from 'vue-router'
// import 'bootstrap/dist/css/bootstrap.css'
// import 'bootstrap-vue/dist/bootstrap-vue.css'

const Home = () => import('@/components/Home')
const Search = () => import('@/components/Search')
const info = () => import('@/components/info')
const info_tissue = () => import('@/components/info_tissue')
const documentation = () => import('@/components/documentation')
const detail = () => import('@/components/detail')
const Statistics = () => import('@/components/Statistics')
const annotation = () => import('@/components/annotation')
const gene = () => import('@/components/gene')

// import Home from '@/components/Home'
// import Search from '@/components/Search'
// import info from '@/components/info'
// import Documentation from '@/components/Documentation'
// import detail from '@/components/detail'
// import Statistics from '@/components/Statistics'


Vue.use(Router)


export default new Router({
    routes: [
        {
            path: '/',
            redirect: '/search',
        },
        {
            path: '/Home',
            name: 'Home',
            component: Home,
            redirect: '/search',
            children: [
                {
                    path: '/search',
                    name: 'Search',
                    component: Search
                },
                {
                    path: '/info',
                    name: 'info',
                    component: info,
                },
                {
                    path: '/info_tissue/:tissue',
                    name: 'info_tissue',
                    component: info_tissue,
                },
                {
                    path: '/detail/:DatasetName',
                    name: 'detail',
                    component: detail
                },
                {
                    path: '/gene',
                    name: 'gene',
                    component: gene
                },
                {
                    path: '/annotation',
                    name: 'annotation',
                    component: annotation
                },
                {
                    path: '/documentation',
                    name: 'documentation',
                    component: documentation
                },
                {
                    path: '/Statistics',
                    name: 'Statistics',
                    component: Statistics,
                }
            ]
        },
    ],
    // mode: 'history'
})

