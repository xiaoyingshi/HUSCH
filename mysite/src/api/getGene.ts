import request from "@/utils/request.ts"

export function searchGene(params) {
  return request({
    url: '/detail_gene',
    method: 'post',
    data: params
  })
}


export function getGene(params) {
    return request({
      url: '/get_gene',
      method: 'post',
      data: params
    })
  }


export function getViolin(params) {
    return request({
      url: '/violin_gene',
      method: 'post',
      data: params
    })
  }


  export function getJson(params) {
    return request({
      url: '/get_json',
      method: 'post',
      data: params
    })
  }

export function getHeatmap(params) {
    return request({
        url: '/get_heatmap',
        method: 'post',
        data: params
    })
}

