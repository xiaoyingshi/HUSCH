<template>
  <div id='stat'>
    
    <iframe width="1200" height="900" src="https://datastudio.google.com/embed/reporting/8990ef42-4376-4eb2-b11c-4104cf050736/page/qiZoC" frameborder="0" style="border:0" allowfullscreen></iframe>
  </div>



</template>


<script>
  export default {
    data () {
    return {
      
    }
  },
  mounted () {
    // 这里就拿到了iframe的对象
    console.log(this.$refs.iframe)
    // 这里就拿到了iframe的window对象
    console.log(this.$refs.iframe.contentWindow)
  }
}

    

</script>


<style lang="scss" scoped>
</style>