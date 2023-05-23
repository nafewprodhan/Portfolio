import scrapy
import logging

class CountriesSpider(scrapy.Spider):
    name = "countries"
    allowed_domains = ["www.worldometers.info"]    # Here don't put slash at the end otherwise error can be occured
    start_urls = ["https://www.worldometers.info/world-population/population-by-country/"]

    def parse(self, response):
        
        # here because we are going to open the link also so we need the whole a element with text and href because we will customize the response we will not use getall() method here
        countries = response.xpath("//td/a")

        for country in countries:
            # the response object in countries will be stored all the selected path now individual selected path also have the xpath method where we select the elements by starting like this .// 
            # so to select the name
            name = country.xpath('.//text()').get()
            link = country.xpath('.//@href').get()
    
            # yield {
            #     "country_name": name,
            #     "country_link": link
            # }


            # But here we are not only crawling one single page but also the corresponding countries content so we need to get request from each link and the content from it to do that:
            # yield scrapy.Request(url=link) # but the problem is the link doesn't contain the domain or type of protocols version but only the path
            #To fix this

            # absolute_url = f"https://www.worldometers.info{link}"
            # But still this isn't the best practice, so other way
            # absolute_url = response.urljoin(link)
            # yield scrapy.Request(absolute_url)

            # Also I can do like this at the same time in one command to successfully request each link
            # yield response.follow(url=link)

            # but only sending the request and getting the response is not enough. we have to capture the request so we can scrape -- where we gonna catch the response? 

            # For that we need to create a method which can catch response
            # as names parsed here we have to send that as meta data
            yield response.follow(url=link, callback=self.parse_country, meta={'country_name':name})
            
    # as scrapy sends the responses to each contry link the responses will be sent to the parse_country() method
    def parse_country(self, response):
        # look for response url
        # logging.info(response.url)
        name = response.request.meta['country_name']
        rows = response.xpath("(//table[@class='table table-striped table-bordered table-hover table-condensed table-list'])[1]/tbody/tr")
        for row in rows:
            year = row.xpath(".//td[1]/text()").get()
            population = row.xpath(".//td[2]/strong/text()").get()
            yearly_changes = row.xpath(".//td[3]/text()").get()
            fertility_rate = row.xpath(".//td[7]/text()").get()

            yield {
                "country_name": name,
                "year": year,
                "population": population,
                "yearly_changes": yearly_changes,
                "fertility_rate": fertility_rate
            }
        

# For data extraction -- command
# scrapy crawl countries -o population_dataset.json
# scrapy crawl countries -o population_dataset.csv
# scrapy crawl countries -o population_dataset.xml           


# clearer version of the template

import scrapy

class CountriesSpider(scrapy.Spider):
    name = "countries"
    allowed_domains = ["www.worldometers.info"]
    start_urls = ["https://www.worldometers.info/world-population/population-by-country/"]

    def parse(self, response):
        
        countries = response.xpath("//td/a")

        for country in countries:
            name = country.xpath('.//text()').get()
            link = country.xpath('.//@href').get()

            yield response.follow(url=link, callback=self.parse_country, meta={'country_name':name})
            
    def parse_country(self, response):
        name = response.request.meta['country_name']
        rows = response.xpath("(//table[@class='table table-striped table-bordered table-hover table-condensed table-list'])[1]/tbody/tr")
        for row in rows:
            year = row.xpath(".//td[1]/text()").get()
            population = row.xpath(".//td[2]/strong/text()").get()
            yearly_changes = row.xpath(".//td[3]/text()").get()
            fertility_rate = row.xpath(".//td[7]/text()").get()

            yield {
                "country_name": name,
                "year": year,
                "population": population,
                "yearly_changes": yearly_changes,
                "fertility_rate": fertility_rate
            }