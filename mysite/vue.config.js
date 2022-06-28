const bootstrapSassAbstractsImports = require('vue-cli-plugin-bootstrap-vue/sassAbstractsImports.js')

const webpack = require('webpack')
const CompressionWebpackPlugin = require('compression-webpack-plugin')
const productionGzipExtensions = ['js', 'css']
const isProduction = process.env.NODE_ENV === 'production'
const BundleAnalyzerPlugin = require('webpack-bundle-analyzer').BundleAnalyzerPlugin


module.exports = {
    
    // 输出目录
    // plugins: [
    //     new webpack.IgnorePlugin(/^\.\/locale$/, /moment$/),

    //     // 下面是下载的插件的配置
    //     new CompressionWebpackPlugin({
    //         algorithm: 'gzip',
    //         test: new RegExp('\\.(' + productionGzipExtensions.join('|') + ')$'),
    //         threshold: 10240,
    //         minRatio: 0.8
    //     }),
    //     new webpack.optimize.LimitChunkCountPlugin({
    //         maxChunks: 5,
    //         minChunkSize: 100
    //     }),
    //     new BundleAnalyzerPlugin()
    // ],

    // chainWebpack: (config) => {
    //     // 只在生产环境使用 cdn
    //     if (process.env.NODE_ENV === "production") {
    //         // 忽略 vue 和 moment 这两个模块
    //         config.externals({
    //             // 'vue': "Vue",
    //             'axios': "axios",
    //             'element-ui': 'ELEMENT'
    //         });

    //         // 修改 HtmlWebpackPlugin 插件参数，植入 cdns 这个模板参数，值为 Vue3 和 Moment.js 的 cdn 链接
    //         config.plugin("html").tap((args) => {
    //             args[0].cdns = `
    //             // <script src="https://cdn.jsdelivr.net/npm/vue@2.6.14/dist/vue.esm.browser.js"></script>
    //             <script src="https://unpkg.com/vue-router@3.4.3/dist/vue-router.js"></script>
    //             <script src="https://cdn.staticfile.org/axios/0.21.4/axios.min.js"></script>
    //             <script src="https://unpkg.com/element-ui@2.13.0/lib/index.js"></script>
    //                                 `;
    //             return args;
    //         });
    //     }
    // },

    assetsDir: 'static',
    publicPath: '/',
    css: {
        loaderOptions: {
            sass: {
                additionalData: bootstrapSassAbstractsImports.join('\n')
            },
            scss: {
                additionalData: [...bootstrapSassAbstractsImports, ''].join(';\n')
            }
        }
    },
    configureWebpack: {
        plugins: [
            // new BundleAnalyzerPlugin({
            //     analyzerPort: 4000,
            // })   // 打包使用分析插件
            // ,
            new webpack.NormalModuleReplacementPlugin(/element-ui\/lib\/locale\/lang\/zh-CN/, 'element-ui/lib/locale/lang/en')

        ],
        resolve: {
            alias: {
                'assets': '@/assets',
                // 'common': '@/common',
                'components': '@/components',
                'styles': '@/styles',
                'api': '@/api',

                // 'views': '@/views',
            }
        },

        performance: {
            hints: 'warning',
            //入口起点的最大体积
            maxEntrypointSize: 50000000,
            //生成文件的最大体积
            maxAssetSize: 30000000,
            //只给出 js 文件的性能提示
            assetFilter: function (assetFilename) {
                return assetFilename.endsWith('.js');
            }
        }
    },
    devServer: {
        port: 8080,
        host: "localhost",
        https: false,
        open: false,
        proxy: {
            '/api': {
                target: 'http://127.0.0.1:8000/api/', //接口域名
                // target: 'http://39.101.160.221:8000/api/', //接口域名

                changeOrigin: true, //是否跨域
                pathRewrite: {
                    '^/api': '',//需要rewrite重写的
                }
            }
        }
    },
    //webpack配置
}


