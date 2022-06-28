// "use strict";

// import axios, { AxiosRequestConfig, AxiosResponse, AxiosInstance } from "axios";
// import { Message } from 'element-ui';

// const config = {
//   baseURL: process.env.VUE_APP_BASE_API,
//   timeout: 360 * 1000, // Timeout
//   // withCredentials: true, // Check cross-site Access-Control
// };

// const axiosInstance: AxiosInstance = axios.create(config);

// // Add request interceptor
// axiosInstance.interceptors.request.use(
//   (config: AxiosRequestConfig) => {
//     return config;
//   },
//   (error: any) => {
//     console.log(error);
//     Promise.reject(error);
//   }
// );

// // Add a response interceptor
// axiosInstance.interceptors.response.use(
//   (response: AxiosResponse) => {
//     const res = response.data
//     if (res.code !== 200) {
//       Message({
//         message: res.message,
//         type: 'error',
//         duration: 3 * 1000
//       })
//     } else {
//       return res
//     }
//   },
//   (error: any) => {
//     console.log(error);
//     Message({
//       message: error.message,
//       type: 'error',
//       duration: 3 * 1000
//     })
//     return Promise.reject(error)
//   }
// );

// export default axiosInstance;
