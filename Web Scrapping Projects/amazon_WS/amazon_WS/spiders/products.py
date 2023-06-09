import scrapy


class ProductsSpider(scrapy.Spider):
    name = "products"
    allowed_domains = ["www.ebay.com"]
    # start_urls = ["https://www.ebay.com/b/Samsung/bn_21834655"]     # sending request to the targetted web

    def start_requests(self):
        yield scrapy.Request(url="https://www.ebay.com/b/Samsung/bn_21834655", callback=self.parse, headers={
            'User-Agent': "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/111.0.0.0 Safari/537.36"
        })

    # after sending the request we get a response from the web
    def parse(self, response):
        # we store all the targetted objects
        products = response.xpath('//ul[@class="b-list__items_nofooter srp-results srp-grid"]/li')
        # to yield or extract the targetted elements in each objects
        for product in products:
            yield {
                'title' : product.xpath('.//div/div[@class="s-item__info clearfix"]/a/h3/text()').get(),
                'url' : response.urljoin(product.xpath('.//div/div[@class="s-item__info clearfix"]/a/@href').get()),
                # 'price' : product.xpath('.//').get()
                # 'sold' : product.xpath('.//').get()
                'User-Agent': response.request.headers['User-Agent']  #just to check
            }

        # If there are any pagination or next page link available after yielding all the targetted objects in the page
        next_page = response.xpath('//div[@class="b-pagination"]/nav/a/@href').get()

        if next_page:
            # If available, send a request to the next page url and call the parse method to do the same task until there are no pagination button: so the new request will catch new response for new page and so on.
            yield scrapy.Request(url=next_page, callback=self.parse, headers={
                    'User-Agent': "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/111.0.0.0 Safari/537.36"
            })